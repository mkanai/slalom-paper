import argparse
import hail as hl
from fm_insights.utils import checkpoint_tmp, register_log


pip_bin_breaks = [0, 0.01, 0.1, 0.5, 0.9, 1.0]
q_bin_breaks = [0, 1e-4, 5e-4, 1e-3, 5e-3, 0.01, 1]


def agg_bin(expr, breaks=pip_bin_breaks, expr_group=None, f=lambda x: hl.agg.count(), include_highest=True):
    if expr_group is None:
        expr_group = expr
    if expr_group.dtype == hl.tfloat64:
        texpr_group = hl.tfloat64
    elif expr_group.dtype == hl.tint64:
        texpr_group = hl.tint64
    include_highest = hl.literal(include_highest)

    _bin_idx_f = hl.experimental.define_function(
        lambda v: (
            hl.case()
            .when(v < 0, -1)
            .when((breaks[0] <= v) & (v < breaks[1]), 0)
            .when((breaks[1] <= v) & (v < breaks[2]), 1)
            .when((breaks[2] <= v) & (v < breaks[3]), 2)
            .when((breaks[3] <= v) & (v < breaks[4]), 3)
            .when(
                include_highest,
                (
                    hl.case()
                    .when((breaks[4] <= v) & (v <= breaks[5]), 4)
                    .when(v > breaks[5], 5)
                    .default(hl.missing(hl.tint32))
                ),
            )
            .when(
                ~include_highest,
                (
                    hl.case()
                    .when((breaks[4] <= v) & (v < breaks[5]), 4)
                    .when(v >= breaks[5], 5)
                    .default(hl.missing(hl.tint32))
                ),
            )
            .default(hl.missing(hl.tint32))
        ),
        texpr_group,
    )

    ret = hl.agg.group_by(_bin_idx_f(expr_group), f(expr))
    return hl.range(0, 5).map(lambda i: ret.get(i))


def main(args):
    suffix = "_grch38" if args.grch38 else ""
    reference_genome = "GRCh38" if args.grch38 else "GRCh37"
    chr_prefix = "chr" if args.grch38 else ""
    basedir = "gs://meta-finemapping-simulation/new_simulations"

    if args.munge:
        ht = hl.import_table(
            f"{basedir}/abf{suffix}/*.ABF.snp.bgz",
            impute=True,
            types={"t_dentist": hl.tfloat64, "nlog10p_dentist": hl.tfloat64},
        )
        ht = ht.annotate(
            pheno=hl.int32(ht.trait.split("\\.")[1][5:]),
            config=hl.int32(ht.trait.split("\\.")[2][6:]),
            **hl.parse_variant(ht.variant, reference_genome=reference_genome),
        )
        ht = ht.key_by("pheno", "config", "locus", "alleles")
        ht = checkpoint_tmp(ht)

        # get unique regions
        ht_region = ht.group_by("pheno", "config", "region").aggregate(
            gamma=hl.agg.any(ht.gamma),
            gamma_variant=hl.agg.take(ht.variant, 1, ordering=-ht.gamma)[0],
            gamma_prob=hl.agg.take(ht.prob, 1, ordering=-ht.gamma)[0],
            max_prob=hl.agg.max(ht.prob),
            max_nlog10p_dentist=hl.agg.max(ht.nlog10p_dentist),
            lead_z=hl.agg.filter(ht.lead_variant, hl.agg.collect(ht.beta / ht.se))[0],
            cs_variant=hl.agg.filter(ht.cs, hl.agg.collect_as_set(ht.variant)),
            cs_99_variant=hl.agg.filter(ht.cs_99, hl.agg.collect_as_set(ht.variant)),
            variant_prob=hl.agg.collect_as_set(hl.struct(variant=ht.variant, prob=ht.prob)),
        )

        ht_region = ht_region.annotate(
            gamma_variant=hl.or_missing(ht_region.gamma, ht_region.gamma_variant),
            gamma_prob=hl.or_missing(ht_region.gamma, ht_region.gamma_prob),
            max_prob_variants=ht_region.variant_prob.filter(lambda x: x.prob == ht_region.max_prob).variant,
        )
        ht_region = ht_region.annotate(
            abf_success=(ht_region.max_prob_variants.contains(ht_region.gamma_variant)) & ht_region.gamma,
            abf_success_cs=ht_region.cs_variant.contains(ht_region.gamma_variant),
            abf_success_cs_99=ht_region.cs_99_variant.contains(ht_region.gamma_variant),
            interval=hl.parse_locus_interval(
                chr_prefix + ht_region.region[3:].replace(":0-", ":1-"), reference_genome=reference_genome
            ),
        )
        ht_region = ht_region.drop("variant_prob")
        ht_region = ht_region.key_by("pheno")
        ht_region = ht_region.checkpoint(f"{basedir}/sim{suffix}.ALL.region.ht", overwrite=args.overwrite)

        # get true causal variants
        liftover = ".liftover" if args.grch38 else ""
        ht_beta = hl.import_table(
            f"{basedir}/true_pheno/pheno*.beta{liftover}.txt", impute=True, source_file_field="source",
        )
        ht_beta = ht_beta.select(
            pheno=hl.int32(ht_beta.source.split("/")[-1].split("\\.")[0][5:]),
            **hl.parse_variant(ht_beta.variant, reference_genome=reference_genome),
        )
        ht_beta = ht_beta.key_by("pheno")
        ht_beta = ht_beta.join(ht_region.key_by("pheno"), "inner")
        ht_beta = ht_beta.filter(ht_beta.interval.contains(ht_beta.locus))
        ht_beta = ht_beta.key_by("pheno", "config", "locus", "alleles").select(gamma_from_region=True)
        ht_beta = checkpoint_tmp(ht_beta)

        # add missing true causal variants
        ht = ht.join(ht_beta, "outer")
        ht = ht.annotate(
            pheno_group=(ht.pheno - 1) // 100 + 1,
            config_group=hl.case()
            .when((13 <= ht.config) & (ht.config <= 17), 1317)  # EUR array mixed, 1000GP3
            .when((18 <= ht.config) & (ht.config <= 22), 1822)  # EUR Omni25, imputaion mixed
            .when((23 <= ht.config) & (ht.config <= 27), 2327)  # EUR array mixed, imputaion mixed
            .when((28 <= ht.config) & (ht.config <= 32), 2832)  # multi-ancestry array mixed, 1000GP3
            .when((33 <= ht.config) & (ht.config <= 37), 3337)  # multi-ancestry Omni25, imputation mixed
            .when((38 <= ht.config) & (ht.config <= 42), 3842)  # multi-ancestry array mixed, imputation mixed
            .default(ht.config),
        )
        ht = ht.annotate(id=hl.delimit([ht.pheno_group, ht.config_group]))
        ht = ht.order_by("id", "prob").add_index()
        ht = checkpoint_tmp(ht)

        # compute index
        ht_idx = ht.group_by("id").aggregate(
            max_idx=hl.agg.max(hl.or_missing(hl.is_defined(ht.prob), ht.idx)),  # NAs located last
            n=hl.agg.count_where(hl.is_defined(ht.prob)),
        )

        d_max_idx = hl.literal(ht_idx.select("max_idx").to_pandas().set_index("id").T.to_dict("records")[0])
        d_n = hl.literal(ht_idx.select("n").to_pandas().set_index("id").T.to_dict("records")[0])

        # add index
        ht = ht.annotate(rank=d_max_idx[ht.id] - ht.idx)
        ht = ht.annotate(q=hl.float64(ht.rank / d_n[ht.id]))
        ht = ht.checkpoint(f"{basedir}/sim{suffix}.ALL.ht", overwrite=args.overwrite)

    if args.export_metrics:
        ht = hl.read_table(f"{basedir}/sim{suffix}.ALL.ht")
        ht_hist = ht.group_by("pheno_group", "config_group").aggregate(
            prob_bins=agg_bin(ht.prob),
            prob_true_bins=agg_bin(hl.or_missing(ht.gamma, ht.prob)),
            mean_prob=agg_bin(ht.prob, f=hl.agg.stats),
            q_true_bins=agg_bin(hl.or_missing(ht.gamma, ht.q), breaks=q_bin_breaks),
            Ngamma_region=hl.agg.sum(ht.gamma_from_region),
        )
        ht_hist = ht_hist.annotate(
            calibration=ht_hist.prob_true_bins / ht_hist.prob_bins,
            mean_prob=ht_hist.mean_prob.mean,
            mean_prob_sd=ht_hist.mean_prob.stdev,
        )
        ht_hist.export(f"{basedir}/sim{suffix}.ALL.metrics.tsv.bgz")

    if args.export_dentist:
        ht = hl.read_table(f"{basedir}/sim{suffix}.ALL.ht")
        ht_region = hl.read_table(f"{basedir}/sim{suffix}.ALL.region.ht")
        ht_region = ht_region.filter(ht_region.gamma)
        ht_region = ht_region.key_by("pheno", "config", "region")
        ht = ht.key_by("pheno", "config", "region").join(ht_region, "inner")

        ht2 = ht.filter(ht.lead_variant | ht.gamma)
        ht2.select("n_samples", "lead_variant", "gamma", "prob", "p_het").export(
            f"{basedir}/sim{suffix}.ALL.n_samples.tsv.bgz"
        )

        ht2 = ht.select(
            "gamma", "gamma_from_region", "prob", "cs", "cs_99", "lead_variant", "t_dentist", "nlog10p_dentist"
        )
        ht2.filter(ht2.prob > 0.01).export(f"{basedir}/sim{suffix}.ALL.prob001.tsv.bgz")

        _bin_idx_f = hl.experimental.define_function(
            lambda s, e, nbins, binsize, v: (
                hl.case()
                .when(v < s, -1)
                .when(v == s, 0)
                .when(v > e, nbins)
                .default(hl.int32(hl.ceil((v - s) / binsize) - 1))
            ),
            hl.tfloat64,
            hl.tfloat64,
            hl.tint32,
            hl.tfloat64,
            hl.tfloat64,
        )

        ht = ht.filter(~ht.lead_variant)
        ht = ht.annotate(
            r2_bin=_bin_idx_f(0.0, 1.0, 10, 0.1, ht.r ** 2),
            max_prob_bin=_bin_idx_f(0.0, 1.0, 10, 0.1, ht.max_prob),
            nlog10p_dentist_bin=_bin_idx_f(0.0, 15.0, 30, 0.5, ht.nlog10p_dentist),
        )
        ht = ht.group_by("pheno", "config", "region", "r2_bin", "max_prob_bin", "nlog10p_dentist_bin").aggregate(
            n=hl.agg.count(),
            n_above=hl.agg.count_where((ht.r * ht.lead_z - ht.beta / ht.se) > 0),
            n_below=hl.agg.count_where((ht.r * ht.lead_z - ht.beta / ht.se) < 0),
            weighted_n_by_n_samples=hl.agg.sum(ht.n_samples),
            weighted_n_by_dentist=hl.agg.sum(hl.if_else(ht.nlog10p_dentist > 100, 100, ht.nlog10p_dentist)),
            max_nlog10p_dentist=hl.agg.max(ht.nlog10p_dentist),
            max_n_samples=hl.agg.max(ht.n_samples),
            min_n_samples=hl.agg.min(ht.n_samples),
            abf_success=hl.agg.any(ht.abf_success),
            abf_success_cs=hl.agg.any(ht.abf_success_cs),
            abf_success_cs_99=hl.agg.any(ht.abf_success_cs_99),
        )
        ht.export(f"{basedir}/sim{suffix}.ALL.dentist.frac.tsv.bgz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--munge", action="store_true")
    parser.add_argument("--export-metrics", action="store_true")
    parser.add_argument("--export-dentist", action="store_true")
    parser.add_argument("--grch38", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
