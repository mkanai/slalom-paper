import hail as hl
from fm_insights import get_merged_results_path, get_merged_eqtl_path, liftover

ht = hl.read_table(get_merged_results_path("fm_only.max_pip_per_var"))
ht = liftover(ht)
ht = ht.key_by()
ht = ht.select(
    variant=hl.variant_str(ht.locus, ht.alleles),
    variant_b37=hl.variant_str(ht.original_locus, ht.original_alleles),
    max_pip=ht.max_pip,
)
ht.export("gs://meta-finemapping-simulation/fm_only.max_pip_per_var.b38.tsv.bgz")

ht = hl.read_table(get_merged_eqtl_path("merged.max_pip_per_var"))
ht = liftover(ht)
ht = ht.key_by()
ht = ht.select(
    variant=hl.variant_str(ht.locus, ht.alleles),
    variant_b37=hl.variant_str(ht.original_locus, ht.original_alleles),
    max_pip=ht.max_pip,
)
ht.export("gs://meta-finemapping-simulation/merged.max_pip_per_var.b38.tsv.bgz")
