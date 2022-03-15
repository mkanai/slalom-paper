import argparse
import hail as hl
from fm_insights import checkpoint_tmp, get_merged_results_path, liftover, register_log, POPS

bucket = "meta-finemapping-simulation/biobank_finemapping"


def get_raw_susie_path(pop, trait="*", sim=False):
    suffix = "_prs_pheno" if sim else ""
    raw_susie_path = {
        "BBJ": f"gs://{bucket}/bbj/BBJ.{trait}.SuSiE.snp.bgz",
        "UKBB": f"gs://{bucket}/ukbb/UKBB.{trait}.SuSiE.snp.bgz",
        "FG": f"gs://{bucket}/fg_r6/{trait}.SUSIE.snp.bgz",
    }
    return raw_susie_path[pop]


def get_raw_finemap_path(pop, trait="*", sim=False):
    suffix = "_prs_pheno" if sim else ""
    raw_finemap_path = {
        "BBJ": f"gs://{bucket}/bbj/BBJ.{trait}.FINEMAP.snp.bgz",
        "UKBB": f"gs://{bucket}/ukbb/UKBB.{trait}.FINEMAP.snp.bgz",
        "FG": f"gs://{bucket}/fg_r6/{trait}.FINEMAP.snp.bgz",
    }
    return raw_finemap_path[pop]


def bbj_filter_f(ht):
    # remove imputed SVs for now (e.g., <INS:ME:ALU>)
    return ~ht.variant.contains("<")


def fg_select_f(ht):
    expr = {
        # "cohort": ht.cohort,
        "trait": ht.trait,
        "region": ht.region,
        "variant": "chr" + ht.v,
        # "variant_b37": ht.variant_b37,
        "prob": ht.prob,
        "cs": ht.cs,
        "mean": ht.mean,
        "sd": ht.sd,
        # "need_to_flip_beta": ht.need_to_flip_beta,
    }
    if "low_purity" in ht.row:
        expr["low_purity"] = ht.low_purity
    if "alpha1" in ht.row:
        expr = {**expr, **{f"alpha{i}": ht[f"alpha{i}"] for i in range(1, 11)}}
    return expr


def read_snp(path, method, trait_dict, filter_f=None, select_f=None, reference_genome="GRCh37"):
    ht = hl.import_table(path, min_partitions=200, impute=True)
    if filter_f is not None:
        ht = ht.filter(filter_f(ht))
    if select_f is not None:
        ht = ht.select(**select_f(ht))

    ht = ht.annotate(**hl.parse_variant(ht.variant, reference_genome=reference_genome))

    # map trait_cohort to trait
    ht = ht.annotate(trait=trait_dict[ht.trait], trait_cohort=ht.trait)
    ht = ht.key_by("locus", "alleles", "trait")

    # remove variants in MHC region
    chrom6 = "6" if reference_genome == "GRCh37" else "chr6"
    MHC_region = hl.interval(
        hl.locus(chrom6, 25000000, reference_genome=reference_genome),
        hl.locus(chrom6, 36000000, reference_genome=reference_genome),
        includes_end=True,
    )
    ht = ht.annotate(in_MHC=MHC_region.contains(ht.locus))
    ht = ht.filter(~ht.in_MHC)

    if "low_purity" in ht.row:
        ht = ht.annotate(cs=hl.if_else((ht.cs == -1) | (ht.low_purity == 1), -1, ht.cs))

    if "alpha1" in ht.row:
        susie_alpha_expr = {"alpha": hl.array([ht[f"alpha{i}"] for i in range(1, 11)])}
    else:
        susie_alpha_expr = {"alpha": hl.missing(hl.tarray(hl.tfloat64))}

    if "lbf_variable1" in ht.row:
        susie_lbf_variable_expr = {"lbf_variable": hl.array([ht[f"lbf_variable{i}"] for i in range(1, 11)])}
    else:
        susie_lbf_variable_expr = {"lbf_variable": hl.missing(hl.tarray(hl.tfloat64))}

    ht = ht.select(
        method=method,
        region=ht.region,
        pip=ht.prob,
        cs_id=ht.cs,
        beta_posterior=ht.mean,
        sd_posterior=ht.sd,
        **susie_alpha_expr,
        **susie_lbf_variable_expr,
    )
    return ht


def munge_results(pop, trait="*", filter_f=None, select_f=None, sim=False, overwrite=False):
    trait_dict = hl.literal(
        {
            "BBJ": {"CHF_Combined": "HF", "Hyperuricemia": "Gout"},
            "FG": {"C3_THYROID_GLAND": "ThC", "GOUT": "Gout", "I9_HEARTFAIL_NS": "HF", "I9_VTE": "VTE"},
            "UKBB": {"ThC": "ThC", "Gout": "Gout", "HF": "HF", "VTE": "VTE"},
        }[pop]
    )
    reference_genome = "GRCh37" if pop != "FG" else "GRCh38"

    ht_susie = read_snp(
        get_raw_susie_path(pop, trait, sim=sim),
        "SUSIE",
        trait_dict,
        filter_f=filter_f,
        select_f=select_f,
        reference_genome=reference_genome,
    )
    ht_susie = checkpoint_tmp(ht_susie)

    ht_finemap = read_snp(
        get_raw_finemap_path(pop, trait, sim=sim),
        "FINEMAP",
        trait_dict,
        filter_f=filter_f,
        select_f=select_f,
        reference_genome=reference_genome,
    )
    ht_finemap = checkpoint_tmp(ht_finemap)

    ht = ht_susie.union(ht_finemap)
    ht = checkpoint_tmp(ht)

    ht = ht.collect_by_key()
    ht = ht.transmute(
        method="Average",
        region=ht.values.region[0],
        pip=(
            hl.case()
            .when(hl.len(ht.values) == 1, ht.values.pip[0])
            .when(hl.is_missing(ht.values.pip[0]), ht.values.pip[1])
            .when(hl.is_missing(ht.values.pip[1]), ht.values.pip[0])
            .when(hl.abs(ht.values.pip[0] - ht.values.pip[1]) < 0.05, hl.mean(ht.values.pip))
            .or_missing()
        ),
        finemap=hl.struct(
            pip=ht.values.pip[ht.values.method.index("FINEMAP")],
            cs_id=ht.values.cs_id[ht.values.method.index("FINEMAP")],
            beta_posterior=ht.values.beta_posterior[ht.values.method.index("FINEMAP")],
            sd_posterior=ht.values.sd_posterior[ht.values.method.index("FINEMAP")],
        ),
        susie=hl.struct(
            pip=ht.values.pip[ht.values.method.index("SUSIE")],
            cs_id=ht.values.cs_id[ht.values.method.index("SUSIE")],
            beta_posterior=ht.values.beta_posterior[ht.values.method.index("SUSIE")],
            sd_posterior=ht.values.sd_posterior[ht.values.method.index("SUSIE")],
            alpha=ht.values.alpha[ht.values.method.index("SUSIE")],
            lbf_variable=ht.values.lbf_variable[ht.values.method.index("SUSIE")],
        ),
    )
    ht.describe()
    ht = checkpoint_tmp(ht)

    ht = liftover(ht)
    if pop != "FG":
        ht = ht.key_by(
            "trait",
            variant=hl.variant_str(ht.locus, ht.alleles),
            variant_b37=hl.variant_str(ht.original_locus, ht.original_alleles),
        )
    else:
        ht = ht.key_by(
            "trait",
            variant=hl.variant_str(ht.original_locus, ht.original_alleles),
            variant_b37=hl.variant_str(ht.locus, ht.alleles),
        )

    ht = ht.select(cohort=pop, pip=ht.pip, region=ht.region, cs_id=ht.susie.cs_id)
    ht = ht.checkpoint(f"gs://{bucket}/results/{pop}.ht", overwrite=overwrite)
    return ht


def get_previous_results():
    traits = hl.set(["Asthma", "COPD", "IPF", "Glaucoma", "IS"])

    ht = hl.read_table(get_merged_results_path("fm_only"))
    ht = ht.filter(traits.contains(ht.trait))

    ht = liftover(ht)
    ht = ht.key_by(
        "trait",
        variant=hl.variant_str(ht.locus, ht.alleles),
        variant_b37=hl.variant_str(ht.original_locus, ht.original_alleles),
    )
    ht = ht.select(cohort=ht.cohort, pip=ht.pip, region=ht.region, cs_id=ht.susie.cs_id)
    return ht


def export_max_pip(ht, overwrite=False):
    ht = ht.collect_by_key()
    ht = ht.select(
        max_pip=hl.max(ht.values.pip),
        any_cs=hl.any(hl.map(lambda x: x > 0, ht.values.cs_id)),
        **{
            pop: hl.rbind(
                ht.values.cohort.index(pop), lambda idx: hl.or_missing(hl.is_defined(idx), ht.values.pip[idx]),
            )
            for pop in POPS
        },
    )
    ht.describe()
    ht = ht.checkpoint("gs://meta-finemapping-simulation/bbj_fg_ukbb.shared_trait.max_pip.ht", overwrite=overwrite)
    ht.export("gs://meta-finemapping-simulation/bbj_fg_ukbb.shared_trait.pip.tsv.bgz")
    return ht


def compute_cs_size(ht):
    ht2 = ht.filter(ht.cs_id > 0)
    ht2 = ht2.group_by(ht2.trait, ht2.cohort, ht2.region, ht2.cs_id).aggregate(cs_size=hl.agg.count())
    ht2 = ht2.checkpoint("gs://meta-finemapping-simulation/bbj_fg_ukbb.shared_trait.cs_size.ht")
    return ht2


def main(args):
    if args.munge:
        ht_bbj = munge_results("BBJ", trait="*", filter_f=bbj_filter_f, overwrite=args.overwrite)
        ht_ukbb = munge_results("UKBB", trait="*", overwrite=args.overwrite)
        ht_fg = munge_results("FG", trait="*", select_f=fg_select_f, overwrite=args.overwrite)

        ht = hl.Table.union(*[ht_bbj, ht_ukbb, ht_fg, get_previous_results()])
        ht = ht.checkpoint("gs://meta-finemapping-simulation/bbj_fg_ukbb.shared_trait.pip.ht", overwrite=args.overwrite)
    else:
        ht = hl.read_table("gs://meta-finemapping-simulation/bbj_fg_ukbb.shared_trait.pip.ht")

    if args.export_max_pip:
        export_max_pip(ht, overwrite=args.overwrite)
    if args.export_cs_size:
        ht2 = compute_cs_size(ht)

        v = hl.import_table("gs://meta-finemapping-simulation/gbmi-all-biobank-meta/lead_variants.txt", impute=True)
        v = v.key_by("variant_b37", "trait")
        ht = ht.key_by("variant_b37", "trait")
        ht = ht.join(v, "inner")
        ht = ht.annotate(cs_size=ht2[ht.trait, ht.cohort, ht.region, ht.cs_id].cs_size)
        ht = ht.key_by()
        ht.select("trait", "cohort", "variant", "cs_size").export(
            "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/lead_cs_size.txt.bgz"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--munge", action="store_true")
    parser.add_argument("--export-max-pip", action="store_true")
    parser.add_argument("--export-cs-size", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
