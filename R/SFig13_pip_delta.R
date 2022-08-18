library(ggplot2)
library(dplyr)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")


df.delta <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi.bbj_fg_ukbb.delta.tsv.bgz")
df.delta2 <-
  dplyr::filter(df.delta, !is.na(prob) & !is.na(max_pip))

plt <-
  dplyr::filter(df.delta, prediction != "NA" & (prob > 0.001 | max_pip > 0.001)) %>%
  dplyr::mutate(prediction = factor(prediction, levels = c("SL", "NSL"))) %>%
  dplyr::group_split(trait) %>%
  purrr::map(function(data) {
    trait <- data$trait[1]
    data <- dplyr::mutate(
      data,
      rsid = dplyr::case_when(
        variant == "chr9:6197392:T:C" ~ "rs1888909",
        variant == "chr5:14610200:C:G" ~ "rs16903574",
        variant == "chr5:132660151:T:C" ~ "rs1295686",
        variant == "chr1:152313385:G:A" ~ "rs61816761",
        is.na(rsid) ~ variant,
        TRUE ~ rsid
      ),
      label = ifelse(prediction == "NSL" &
        abs(delta) > 0.5, rsid, "")
    )

    data.mean <- dplyr::group_by(data, prediction) %>%
      dplyr::summarize(x = mean(delta, na.rm = T))

    ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggrastr::rasterize(geom_point(
        aes(delta, prediction, color = prediction),
        shape = "|",
        data = dplyr::filter(data, abs(delta) <= 0.5)
      ), dpi = 300) +
      geom_point(aes(delta, prediction, color = prediction),
        data = dplyr::filter(data, abs(delta) > 0.5)
      ) +
      geom_point(
        aes(x, prediction),
        color = "black",
        fill = "white",
        size = 2,
        shape = 23,
        data = data.mean
      ) +
      ggrepel::geom_text_repel(
        aes(delta, prediction, label = label),
        size = 2,
        point.size = 4,
        # direction = "x",
        nudge_y = 0.2,
        force = 20,
        segment.size = unit(0.25, "mm"),
        min.segment.length = 0.2,
        segment.curvature = -0.1,
        segment.ncp = 3,
        segment.angle = 20,
        seed = 1000,
        ylim = c(0, Inf),
        data = data
      ) +
      scale_color_manual(values = c(BuenColors::jdb_palette("brewer_celsius")[c(8, 2)], "grey50")) +
      scale_x_continuous(limits = c(-1, 1)) +
      scale_y_discrete(limits = c("SL", "NSL"), expand = expansion(c(0.5, 1))) +
      locusviz::get_default_theme(
        hide.ytitle = TRUE,
        hide.xlab = trait != "VTE",
        legend.position = "none"
      ) +
      labs(
        title = trait,
        x = "ΔPIP (GBMI - individual biobank)"
      ) +
      coord_cartesian(clip = "off")
  }) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig13_pip_delta.pdf",
  plt,
  base_height = 7.2,
  base_width = 7.2,
  device = cairo_pdf
)