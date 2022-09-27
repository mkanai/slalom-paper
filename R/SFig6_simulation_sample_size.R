library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

df.n_samples <-
  dplyr::bind_rows(
    rgsutil::read_gsfile(
      "gs://meta-finemapping-simulation/new_simulations/sim.ALL.n_samples.tsv.bgz"
    ),
    rgsutil::read_gsfile(
      "gs://meta-finemapping-simulation/new_simulations/sim_grch38.ALL.n_samples.tsv.bgz"
    )
  )

df.ratio <-
  dplyr::group_by(df.n_samples, pheno, config, region) %>%
  dplyr::filter(n() > 1 & max(prob) > 0.9 &
    config %in% c(seq(18, 27), seq(33, 42))) %>%
  dplyr::summarize(
    ratio = max(n_samples) / min(n_samples),
    n = n(),
    min_p_het = min(p_het)
  ) %>%
  dplyr::ungroup()

summary(df.ratio)
sum(df.ratio$ratio < 1) / nrow(df.ratio)

plt <- plot_sample_size_ratio_distribution(df.ratio, fill = BuenColors::jdb_palette("brewer_blue")[2])

plt

cowplot::save_plot(
  "figures/SFig6_simulation_sample_size_ratio.pdf",
  plt,
  base_height = 3.6,
  base_width = 3.6,
  device = cairo_pdf
)
