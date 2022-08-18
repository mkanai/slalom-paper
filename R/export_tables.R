library(dplyr)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

############################################################################
# STable 5: GWAS catalog summary
################################################################################
df.final <- read.table(
  "~/src/github.com/mkanai/slalom-paper/metadata/gwas_catalog_full_sumstats.tsv",
  T,
  sep = "\t",
  quote = ""
) %>%
  dplyr::left_join(
    read.table(
      "~/src/github.com/mkanai/slalom-paper/metadata/gwas_catalog_meta.tsv",
      T
    )
  )
succeeded <- read.table("~/src/github.com/mkanai/slalom-paper/metadata/gwas_catalog_succeeded.txt")[, 1]

dplyr::filter(df.final, STUDY.ACCESSION %in% succeeded &
  is_meta_analysis) %>%
  dplyr::select(-(path:is_meta_analysis)) %>%
  dplyr::arrange(STUDY.ACCESSION) %>%
  write.table(
    "./tables/STable5_gwas_catalog_metadata.tsv",
    quote = F,
    row.names = F,
    sep = "\t"
  )

################################################################################
# STable 6: GWAS Catalog SLALOM results
################################################################################
munge_slalom_prediction <- function(data) {
  dplyr::mutate(
    data,
    trait = stringr::str_split_fixed(trait, "_", 2)[, 1],
    prediction = factor(
      dplyr::case_when(
        max_pip < 0.1 ~ "NA",
        n_dentist_outlier > 0 ~ "SL",
        TRUE ~ "NSL"
      ),
      levels = c("SL", "NSL", "NA"),
      labels = c("Suspicious loci", "Non-suspicious loci", "NA"),
    ),
    max_pip_nonsyn = ifelse(is.infinite(max_pip_nonsyn), NA, max_pip_nonsyn)
  ) %>%
    dplyr::select(
      trait,
      region,
      prediction,
      n_total,
      n_r2,
      n_dentist_outlier,
      fraction,
      n_nonsyn,
      max_pip,
      max_pip_nonsyn,
      cs_nonsyn,
      cs_99_nonsyn,
      nonsyn_variants
    )
}

rgsutil::read_gsfile(
  "gs://meta-finemapping-simulation/gwas_catalog/gwas_catalog.max_pip.SLALOM.summary.tsv.bgz"
) %>%
  dplyr::left_join(
    read.table(
      "~/src/github.com/mkanai/slalom-paper/metadata/gwas_catalog_meta.tsv",
      T
    ),
    by = c("trait" = "STUDY.ACCESSION")
  ) %>%
  dplyr::filter(is_meta_analysis) %>%
  munge_slalom_prediction() %>%
  write.table(
    "./tables/STable6_gwas_catalog_slalom.tsv",
    quote = F,
    row.names = F,
    sep = "\t"
  )


############################################################################
# STable 7: GBMI metadata
############################################################################
gbmi_description <- c(
  "AAA" = "Abdominal aortic aneurysm",
  "AcApp" = "Acute appendicitis",
  "Appendectomy" = "Appendectomy",
  "Asthma" = "Asthma",
  "COPD" = "Chronic obstructive pulmonary disease",
  "Gout" = "Gout",
  "HCM" = "Hypertrophic cardiomyopathy",
  "HF" = "Heart failure",
  "IPF" =  "Idiopathic pulmonary fibrosis",
  "POAG" = "Primary open angle glaucoma",
  "Stroke" =  "Stroke",
  "ThC" = "Thyroid cancer",
  "UtC" = "Uterine cancer",
  "VTE" = "Venous thromboembolism"
)

r1_dropped_loci <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi_r1_dropped_loci.txt", header = F) %>%
  dplyr::rename(trait = V1, lead_pip_variant = V2)

df.n_loci <- rgsutil::read_gsfile(
  "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi-all-biobank-meta.SLALOM.summary.bgz"
) %>%
  dplyr::mutate(trait = stringr::str_split_fixed(trait, "_", 2)[, 1]) %>%
  dplyr::anti_join(r1_dropped_loci) %>%
  dplyr::group_by(trait) %>%
  dplyr::summarize(n_loci = n())

df.n_samples <- read.table(
  "./metadata/gbmi_sample_size_by_ancestry.txt",
  T,
  sep = "\t"
) %>%
  dplyr::rename(trait = Endpoint)

dplyr::mutate(
  df.n_samples,
  Ancestry = dplyr::case_when(
    X..Total != "-" ~ "Total",
    Ancestry == "ami" ~ "MID",
    TRUE ~ stringr::str_to_upper(Ancestry)
  )
) %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(biobanks = stringr::str_c(sort(unique(Biobank[Ancestry != "Total"])), collapse = ",")) %>%
  dplyr::group_by(trait, Ancestry) %>%
  dplyr::summarize(
    description = gbmi_description[trait][1],
    biobanks = biobanks[1],
    n_cases = sum(Case),
    n_controls = sum(Control),
    n_samples = n_cases + n_controls
  ) %>%
  dplyr::left_join(df.n_loci) %>%
  tidyr::pivot_wider(
    id_cols = c("trait", "description", "n_loci", "biobanks"),
    names_from = "Ancestry",
    values_from = c("n_cases", "n_controls")
  ) %>%
  dplyr::select(
    trait:biobanks,
    dplyr::ends_with("Total"),
    dplyr::ends_with("AFR"),
    dplyr::ends_with("AMR"),
    dplyr::ends_with("EAS"),
    dplyr::ends_with("FIN"),
    dplyr::ends_with("MID"),
    dplyr::ends_with("NFE"),
    dplyr::ends_with("SAS")
  ) %>%
  write.table(
    "./tables/STable7_gbmi_metadata.tsv",
    quote = F,
    row.names = F,
    sep = "\t"
  )

################################################################################
# STable 8: GBMI SLALOM results
################################################################################
rgsutil::read_gsfile(
  "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi-all-biobank-meta.SLALOM.summary.bgz"
) %>%
  dplyr::mutate(trait_short = stringr::str_split_fixed(trait, "_", 2)[, 1]) %>%
  dplyr::anti_join(r1_dropped_loci, by = c("trait_short" = "trait", "lead_pip_variant")) %>%
  dplyr::select(-trait_short) %>%
  munge_slalom_prediction() %>%
  write.table(
    "./tables/STable8_gbmi_slalom.tsv",
    quote = F,
    row.names = F,
    sep = "\t"
  )