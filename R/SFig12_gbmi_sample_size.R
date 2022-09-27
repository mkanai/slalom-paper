library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

r1_dropped_loci <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi_r1_dropped_loci.txt", header = FALSE) %>%
  dplyr::rename(trait_short = V1, lead_pip_variant = V2)

df.gbmi <- rgsutil::read_gsfile(
  "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi-all-biobank-meta.SLALOM.summary.bgz"
) %>%
  dplyr::mutate(
    prediction = factor(
      dplyr::case_when(
        max_pip < 0.1 ~ "NA",
        n_dentist_outlier > 0 ~ "SL",
        TRUE ~ "NSL"
      ),
      levels = c("SL", "NSL", "NA")
    ),
    trait_short = stringr::str_split_fixed(trait, "_", 2)[, 1]
  ) %>%
  dplyr::anti_join(r1_dropped_loci) %>%
  dplyr::select(-trait_short)

plt <-
  list(
    dplyr::mutate(df.gbmi, ratio = max_neff_r2 / min_neff_r2) %>%
      dplyr::filter(prediction == "SL") %>%
      plot_sample_size_ratio_distribution(
        fill = BuenColors::jdb_palette("brewer_celsius")[8],
        title = "Suspicious loci",
        label.y = 1.05,
        hjust = -0.3,
        log10 = TRUE,
        xlim = c(1, 65)
      ) + labs(x = "Effective sample size ratio (max / min)"),
    dplyr::mutate(df.gbmi, ratio = max_neff_r2 / min_neff_r2) %>%
      dplyr::filter(prediction == "NSL") %>%
      plot_sample_size_ratio_distribution(
        fill = BuenColors::jdb_palette("brewer_celsius")[2],
        title = "Non-suspicious loci",
        label.y = 3,
        log10 = TRUE,
        xlim = c(1, 65)
      ) +
      labs(x = "Effective sample size ratio (max / min)") +
      theme(axis.title.y = element_blank())
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(tag_levels = "a")
plt

cowplot::save_plot(
  "figures/SFig12_gbmi_sample_size_ratio.pdf",
  plt,
  base_height = 3.6,
  base_width = 7.2,
  device = cairo_pdf
)
