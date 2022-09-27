library(ggplot2)
library(dplyr)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

df.delta <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi.bbj_fg_ukbb.delta.tsv.bgz")
df.delta2 <-
  dplyr::filter(df.delta, !is.na(prob) & !is.na(max_pip))

plot_enrichment_delta <- function(data,
                                  title,
                                  threshold = 0.01,
                                  hide.ylab = FALSE,
                                  ylim = c(0.8, 2.1),
                                  colors = unname(annot_colors[c("pLoF", "CRE", "Non-genic")])) {
  pd <- position_dodge(width = 0.75)
  dplyr::mutate(data,
    delta_bin = cut(delta, c(-Inf, -threshold, threshold, 1))
  ) %>%
    dplyr::filter(!is.na(delta) & !is.na(consequence)) %>%
    dplyr::group_by(delta_bin, consequence) %>%
    dplyr::count(delta_bin, consequence) %>%
    dplyr::group_by(delta_bin) %>%
    dplyr::mutate(
      total = sum(n),
      n2 = total - n
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_split(consequence) %>%
    purrr::map_dfr(function(data) {
      print(data)
      pos_bin <- sprintf("(%s,1]", threshold)
      neg_bin <- sprintf("(-Inf,-%s]", threshold)
      null_bin <- sprintf("(-%s,%s]", threshold, threshold)

      m1 <-
        magrittr::extract(data, match(c(null_bin, pos_bin), data$delta_bin), c("n2", "n")) %>%
        as.matrix() %>%
        epitools::riskratio(method = "boot")

      m2 <-
        magrittr::extract(data, match(c(null_bin, neg_bin), data$delta_bin), c("n2", "n")) %>%
        as.matrix() %>%
        epitools::riskratio(method = "boot")

      print(m1)
      print(m2)

      tibble::tibble(
        consequence = data$consequence[1],
        enrichment = c(m1$measure[2, "estimate"], m2$measure[2, "estimate"]),
        lower = c(m1$measure[2, "lower"], m2$measure[2, "lower"]),
        upper = c(m1$measure[2, "upper"], m2$measure[2, "upper"]),
        delta_bin = c(sprintf("ΔPIP > %s", threshold), sprintf("ΔPIP < -%s", threshold))
      )
    }) %>%
    ggplot(
      aes(
        consequence,
        enrichment,
        color = consequence,
        fill = consequence,
        shape = delta_bin
      )
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_point(position = pd) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
      width = 0,
      position = pd
    ) +
    locusviz::get_default_theme(hide.xtitle = TRUE, hide.ylab = hide.ylab) +
    theme(
      legend.margin = margin(),
      legend.key.height = unit(0.3, "cm"),
      plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm")
    ) +
    coord_cartesian(ylim = ylim) +
    labs(
      y = "Enrichment",
      title = title,
      color = "Annotation",
      shape = "ΔPIP"
    ) +
    scale_color_manual(values = colors, guide = NULL) +
    scale_fill_manual(values = colors, guide = NULL) +
    scale_shape_manual(values = c(25, 24)) +
    guides(shape = guide_legend(override.aes = list(fill = "black")))
}

plt <- purrr::imap(c(0.01, 0.05, 0.1), function(threshold, i) {
  dplyr::mutate(
    df.delta2,
    consequence = dplyr::case_when(
      consequence %in% c("pLoF", "Missense", "Synonymous", "UTR5", "UTR3", "Promoter") ~ "Coding|UTR|Promoter",
      is.na(consequence) ~ NA_character_,
      TRUE ~ consequence
    ),
    consequence = factor(
      consequence,
      levels = c("Coding|UTR|Promoter", "CRE", "Non-genic")
    )
  ) %>%
    dplyr::mutate(prob > threshold | max_pip > threshold) %>%
    plot_enrichment_delta(
      title = "All loci",
      threshold = threshold,
      ylim = c(0.8, 2.1),
      hide.ylab = i > 1
    )
}) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig15_pip_delta_enrichment.pdf",
  plt,
  base_height = 2.4,
  base_width = 7.2,
  device = cairo_pdf
)
