library(ggplot2)
library(dplyr)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

df.gbmi <- munge_gbmi_sumstats(gbmi_traits["Asthma"], "chr9:3629917-7270078")

df.afr <-
  dplyr::filter(df.gbmi, gnomad_lead_r_afr**2 > 0.6) %>%
  dplyr::summarize(
    min = min(position),
    max = max(position),
    diff = max - min
  )

df.nfe <- dplyr::filter(df.gbmi, gnomad_lead_r_nfe**2 > 0.6) %>%
  dplyr::summarize(
    min = min(position),
    max = max(position),
    diff = max - min
  )

lead_pos <- 6197392
window <- 50000
start <- lead_pos - window
end <- lead_pos + window

background.layers <- list(
  geom_rect(
    aes(
      xmin = min,
      xmax = max,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = locusviz::get_gnomad_colors()["nfe"],
    alpha = 0.05,
    data = df.nfe
  ),
  geom_rect(
    aes(
      xmin = min,
      xmax = max,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = locusviz::get_gnomad_colors()["afr"],
    alpha = 0.1,
    data = df.afr
  )
)


plt <-
  dplyr::filter(df.gbmi, start < position & position < end) %>%
  dplyr::mutate(cs = TRUE) %>%
  dplyr::select(-gnomad_lead_r_amr, -gnomad_lead_r_eas, -gnomad_lead_r_fin) %>%
  locusviz::preprocess(df.gbmi, pip_col = "prob", cs_id_col = "cs") %>%
  locusviz::plot_r2_panel(background.layers = background.layers) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "grey50") +
  geom_point(aes(x = lead_pos, y = 1), shape = 18, size = 3, data = tibble::tibble()) +
  geom_text(
    aes(x = lead_pos, y = 1, label = "rs1888909"),
    size = 2,
    vjust = -1,
    data = tibble::tibble()
  ) +
  geomtextpath::geom_textsegment(
    aes(
      x = min,
      xend = max,
      y = 0.85,
      yend = 0.85,
      label = sprintf("%d kb", round(df.afr$diff / 1000))
    ),
    color = locusviz::get_gnomad_colors()["afr"],
    vjust = -1,
    size = 2,
    arrow = grid::arrow(length = unit(0.05, "in"), ends = "both", type = "open"),
    data = df.afr
  ) +
  geomtextpath::geom_textsegment(
    aes(
      x = min,
      xend = max,
      y = 0.8,
      yend = 0.8,
      label = sprintf("%d kb", round(df.nfe$diff / 1000))
    ),
    color = locusviz::get_gnomad_colors()["nfe"],
    size = 2,
    arrow = grid::arrow(length = unit(0.05, "in"), ends = "both", type = "open"),
    data = df.nfe
  ) +
  labs(x = "Chromosome 9")

cowplot::save_plot(
  "figures/SFig16_rs1888909_ld.pdf",
  plt,
  base_height = 3.6,
  base_width = 7.2,
  device = cairo_pdf
)