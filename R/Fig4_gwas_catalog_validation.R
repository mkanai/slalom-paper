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


plot_new_calibration_nonsyn <- function(data) {
  ggplot() +
    geomtextpath::geom_texthline(
      yintercept = 0,
      label = "Perfect calibration",
      linetype = "dashed",
      color = "grey50",
      size = 2
    ) +
    geom_ribbon(
      aes(
        threshold,
        # mean_prob - calibration,
        # ymin = mean_prob - calibration.lower,
        # ymax = mean_prob - calibration.upper,
        calibration,
        ymin = calibration.lower,
        ymax = calibration.upper,
        fill = prediction
      ),
      alpha = 0.5,
      data = data
    ) +
    geom_path(aes(threshold, calibration,
      # mean_prob - calibration,
      color = prediction
    ),
    data = data
    ) +
    locusviz::get_default_theme() +
    scale_color_manual(values = BuenColors::jdb_palette("brewer_celsius")[c(8, 2)]) +
    scale_fill_manual(
      values = BuenColors::jdb_palette("brewer_celsius")[c(8, 2)],
      guide = "none"
    ) +
    labs(x = "PIP threshold", y = "Mean PIP - % nonsynonymous lead", color = "Prediction") +
    theme(
      legend.position = c(1, 1),
      legend.justification = c(1, 1)
    ) +
    scale_x_continuous(
      limits = c(0.1, 1),
      breaks = c(0.1, 0.25, 0.5, 0.75, 1),
      expand = expansion(0, 0),
      labels = scales::number_format(drop0trailing = TRUE)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::number_format(drop0trailing = TRUE)
    )
}

p <- purrr::map_dfr(seq(0.1, 1, by = 0.01), function(threshold) {
  dplyr::mutate(df.gwas_catalog,
    gamma = max_pip == max_pip_nonsyn
  ) %>%
    dplyr::filter(n_nonsyn > 0 &
      max_pip > threshold &
      prediction != "NA" &
      !is.na(gamma)) %>%
    dplyr::group_by(prediction) %>%
    dplyr::summarize(
      threshold = threshold,
      mean_prob = mean(max_pip, na.rm = T),
      calibration = mean(gamma, na.rm = T),
      calibration.lower = na_binconf(sum(gamma), n(), "Lower"),
      calibration.upper = na_binconf(sum(gamma), n(), "Upper"),
    )
}) %>%
  dplyr::mutate(
    prediction = forcats::fct_recode(
      prediction,
      "Suspicious loci" = "SL",
      "Non-suspicious loci" = "NSL"
    )
  ) %>%
  plot_new_calibration_nonsyn()

plt <- c(
  plot_validataion_wrap(df.gwas_catalog, "GWAS catalog"),
  list(p)
) %>%
  magrittr::extract(c(1, 4, 2, 3)) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 2) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/Fig4_gwas_catalog_validation.pdf",
  plt,
  base_height = 5,
  base_width = 5,
  device = cairo_pdf
)
