library(dplyr)
library(ggplot2)

df <- rgsutil::read_gsfile(
  "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/AAA_Bothsex_inv_var_meta_GBMI_052021.txt.gz",
  select = c("#CHR", "POS", "inv_var_meta_p")
)

plt <-
  dplyr::rename(df, p = inv_var_meta_p, CHR = `#CHR`) %>%
  dplyr::mutate(idx = dplyr::row_number()) %>%
  dplyr::filter(p < 1e-3) %>%
  ggplot(aes(x = idx, y = -log10(p))) +
  ggrastr::rasterize(geom_point(aes(color = factor(CHR %% 2)), size = 0.05), dpi = 900) +
  scale_color_manual(values = c("grey80", BuenColors::jdb_palette("brewer_blue")[7])) +
  locusviz::get_default_theme(hide.xlab = TRUE, hide.ylab = TRUE) +
  theme(
    legend.position = "none",
    axis.line = element_line(size = 0.1),
    axis.ticks = element_blank(),
    plot.background = element_blank()
  ) +
  scale_x_continuous(expand = expansion(0.05, 0)) +
  scale_y_continuous(trans = locusviz::trans_loglog_p(loglog_p = 10), expand = expansion(0, 0))

cowplot::save_plot(
  "figures/Fig1/Fig1_manhattan.pdf",
  plt,
  base_height = 0.6,
  base_width = 1.4,
  device = cairo_pdf,
  bg = "transparent"
)