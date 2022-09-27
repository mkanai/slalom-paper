library(ggplot2)
library(dplyr)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

df.delta <- rgsutil::read_gsfile(
  "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi.bbj_fg_ukbb.delta.tsv.bgz"
)
df.delta2 <-
  dplyr::filter(df.delta, !is.na(prob) & !is.na(max_pip))

dplyr::filter(df.delta2, cs) %>%
  dplyr::select(trait, variant_b37) %>%
  rgsutil::write_gsfile(
    "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/lead_variants.txt",
    overwrite = FALSE
  )

#############

df.cs_size <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/gbmi-all-biobank-meta/lead_cs_size.txt.bgz") %>%
  dplyr::group_by(trait, variant) %>%
  dplyr::summarize(cs_size = min(cs_size, na.rm = T))

df.min_cs_size <- dplyr::inner_join(df.delta2, df.cs_size) %>%
  dplyr::filter(prediction != "NA") %>%
  dplyr::group_by(trait, prediction, region, cs) %>%
  dplyr::summarize(
    cohort = c("GBMI", "Biobanks"),
    cs_size = c(
      n(),
      cs_size[which(lead_variant)]
    )
  ) %>%
  dplyr::filter(all(is.finite(cs_size)))

dplyr::group_by(df.min_cs_size, cohort, prediction) %>%
  dplyr::summarize(
    mean = mean(cs_size, na.rm = T),
    median = median(cs_size, na.rm = T)
  )

plt <-
  purrr::map(c("SL", "NSL"), function(x) {
    data <- dplyr::filter(df.min_cs_size, prediction == x) %>%
      dplyr::mutate(prediction = factor(
        prediction,
        levels = c("NSL", "SL"),
        labels = c("Non-suspicious loci", "Suspicious loci")
      ))

    ggplot(data, aes(cohort, cs_size)) +
      gghalves::geom_half_violin(
        aes(fill = cohort),
        position = position_nudge(x = -0.1),
        draw_quantiles = c(0.25, 0.5, 0.75),
        size = 0.3
      ) +
      gghalves::geom_half_violin(
        aes(fill = cohort),
        position = position_nudge(x = 0.1),
        draw_quantiles = c(0.25, 0.5, 0.75),
        side = "r",
        size = 0.3
      ) +
      geom_line(aes(group = region), alpha = 0.2, size = 0.2) +
      geom_point(size = 0.1, color = "grey50") +
      scale_y_log10(limits = c(1, 1000)) +
      scale_fill_manual(values = BuenColors::jdb_palette("Darjeeling")[c(3, 5)]) +
      locusviz::get_default_theme(
        legend.position = "none",
        hide.xtitle = TRUE,
        hide.ylab = (x == "NSL")
      ) +
      labs(y = "95% CS size", title = data$prediction[1])
  }) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig16_cs_size.pdf",
  plt,
  base_height = 3.2,
  base_width = 7.2,
  device = cairo_pdf
)
