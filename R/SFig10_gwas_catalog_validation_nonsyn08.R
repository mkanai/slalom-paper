library(ggplot2)
library(dplyr)
library(magrittr)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

df.gwas_catalog <- rgsutil::read_gsfile(
  "gs://meta-finemapping-simulation/gwas_catalog/gwas_catalog.max_pip.nonsyn08.SLALOM.summary.tsv.bgz"
) %>%
  dplyr::mutate(prediction = factor(
    dplyr::case_when(
      max_pip < 0.1 ~ "NA",
      n_dentist_outlier > 0 ~ "SL",
      TRUE ~ "NSL"
    ),
    levels = c("SL", "NSL", "NA")
  )) %>%
  dplyr::left_join(
    read.table(
      "~/src/github.com/mkanai/slalom-paper/metadata/gwas_catalog_meta.tsv",
      T
    ),
    by = c("trait" = "STUDY.ACCESSION")
  ) %>%
  dplyr::filter(is_meta_analysis)

plt <-
  dplyr::mutate(
    df.gwas_catalog,
    n_nonsyn = n_nonsyn08,
    max_pip_nonsyn = max_pip_nonsyn08,
    cs_nonsyn = cs_nonsyn08,
    cs_99_nonsyn = cs_99_nonsyn08
  ) %>%
  plot_validataion_wrap("GWAS catalog") %>%
  magrittr::extract2(1)

cowplot::save_plot(
  "figures/SFig10_gwas_catalog_validation_nonsyn08.pdf",
  plt,
  base_height = 2.5,
  base_width = 2.5,
  device = cairo_pdf
)