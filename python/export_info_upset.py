import argparse
import hail as hl
import pandas as pd
import string
from itertools import product
from fm_insights import checkpoint_tmp, register_log


def merge_info(overwrite=False):
    ht = hl.import_table(
        "gs://meta-finemapping-simulation/imputed_info/*.bgz",
        impute=True,
        source_file_field="filename",
        min_partitions=1000,
    )
    ht_exclude = hl.import_table("gs://meta-finemapping-simulation/imputed_info/1000GP.exclude.txt", no_header=True)
    exclude = ht_exclude.aggregate(hl.agg.collect_as_set(ht_exclude.f0), _localize=False)

    ht_liftover = hl.import_table(
        "gs://meta-finemapping-simulation/imputed_info/TopMed_GRCh38.variants.liftover.txt", impute=True
    )
    ht_liftover = ht_liftover.key_by("original_variant")

    ht = ht.filter(~exclude.contains(ht.SNP))
    ht = ht.annotate(
        ancestry=ht.filename.split("\/")[-1].split("_")[0],
        array=ht.filename.split("\/")[-1].split("_")[3],
        imputation=ht.filename.split("\/")[-1].split("_")[4].split("\.")[0],
    )
    ht = ht.select(
        variant=hl.if_else(ht.imputation == "TopMed", ht_liftover[ht.SNP].variant, ht.SNP),
        variant_b38=hl.or_missing(ht.imputation == "TopMed", ht.SNP),
        freq=ht.ALT_Frq,
        rsq=ht.Rsq,
        ancestry=ht.ancestry,
        array=ht.array,
        imputation=ht.imputation,
        label=hl.delimit([ht.ancestry, ht.array, ht.imputation], "_"),
    )
    ht.describe()
    ht = ht.checkpoint("gs://meta-finemapping-simulation/imputed_info/merged_info.ht", overwrite=overwrite)

    # annotate gnomad
    ht_gnomad = hl.read_table(
        "gs://meta-finemapping-simulation/gnomad/gnomad.genomes.r2.1.1.sites.most_severe.max_maf.ht"
    )
    ht_gnomad = ht_gnomad.select("min_maf", "max_maf")
    ht = ht.filter(~ht.variant.contains("<"))
    ht = ht.key_by(**hl.parse_variant(ht.variant))
    ht = ht.join(ht_gnomad, "inner")
    ht = ht.checkpoint("gs://meta-finemapping-simulation/imputed_info/merged_info.gnomad_maf.ht")

    return ht


def export_upset(ht, d, maf_threshold, rsq_threshold, filter_option=""):
    ht = ht.annotate(level=d[ht.label])
    ht = ht.annotate_globals(level_dict=d)
    if filter_option == "":
        ht = ht.filter(
            hl.is_defined(ht.variant) & ((0.5 - hl.abs(0.5 - ht.freq)) > maf_threshold) & (ht.rsq > rsq_threshold)
        )
    elif filter_option == "gnomad_min":
        ht = ht.filter(hl.is_defined(ht.variant) & (ht.min_maf > maf_threshold) & (ht.rsq > rsq_threshold))
    elif filter_option == "gnomad_max_min":
        ht = ht.filter(
            hl.is_defined(ht.variant)
            & (ht.max_maf > maf_threshold)
            & (ht.min_maf <= maf_threshold)
            & (ht.rsq > rsq_threshold)
        )
    maf_prefix = filter_option + "_" if filter_option != "" else ""

    ht = ht.group_by("variant").aggregate(collection=hl.delimit(hl.sorted(hl.agg.collect(ht.level)), ""))
    ht.describe()
    ht = checkpoint_tmp(ht)
    ht = ht.group_by("collection").aggregate(n=hl.agg.count())
    ht.describe()
    ht = checkpoint_tmp(ht)

    maf_str = str(maf_threshold).replace(".", "")
    rsq_str = str(rsq_threshold).replace(".", "")
    ht.export(f"gs://meta-finemapping-simulation/imputed_info/upset.{maf_prefix}maf{maf_str}.rsq{rsq_str}.tsv.bgz")


def main(args):
    if args.merge_info:
        ht = merge_info(overwrite=args.overwrite)
    else:
        ht = hl.read_table("gs://meta-finemapping-simulation/imputed_info/merged_info.gnomad_maf.ht")

    d = hl.dict({x: string.ascii_lowercase[i] for i, x in enumerate(ht.aggregate(hl.agg.collect_as_set(ht.label)))})
    if args.export_level:
        df = pd.DataFrame.from_dict(d.collect()[0], orient="index", columns=["level"])
        with hl.hadoop_open(f"gs://meta-finemapping-simulation/imputed_info/level.txt", "w") as f:
            df.to_csv(f, sep="\t", index_label="label")

    maf_thresholds = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2]
    rsq_thresholds = [0.0, 0.2, 0.4, 0.6, 0.8]
    for maf_threshold, rsq_threshold in product(maf_thresholds, rsq_thresholds):
        export_upset(ht, d, maf_threshold, rsq_threshold, filter_option=args.filter_option)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--merge-info", action="store_true")
    parser.add_argument("--export-upset", action="store_true")
    parser.add_argument("--filter-option", type=str)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
