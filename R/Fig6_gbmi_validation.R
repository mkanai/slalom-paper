library(ggplot2)
library(dplyr)
library(magrittr)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

df.gbmi <- rgsutil::read_gsfile(
  "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi-all-biobank-meta.SLALOM.summary.bgz"
) %>%
  dplyr::mutate(prediction = factor(
    dplyr::case_when(
      max_pip < 0.1 ~ "NA",
      n_dentist_outlier > 0 ~ "SL",
      TRUE ~ "NSL"
    ),
    levels = c("SL", "NSL", "NA")
  ))

plts_validation <- plot_validataion_wrap(df.gbmi, "GBMI")
plts_validation <- purrr::imap(plts_validation, function(x, i) {
  x + labs(tag = letters[i])
})

#########################################

dplyr::mutate(df.gbmi,
  count_lead = max_pip == max_pip_nonsyn
) %>%
  dplyr::filter(n_nonsyn > 0) %$%
  table(prediction, count_lead, useNA = "ifany") %>%
  as.matrix() %>%
  epitools::riskratio()


dplyr::mutate(
  df.gbmi,
  count_lead = lead_max_pip_gwas > 0.9,
  count_cs = cs_max_pip_gwas > 0.9,
  count_cs_99 = cs_99_max_pip_gwas > 0.9
) %$%
  table(prediction, count_lead, useNA = "ifany") %>%
  as.matrix() %>%
  epitools::riskratio()

dplyr::mutate(
  df.gbmi,
  count_lead = lead_max_pip_eqtl > 0.9,
  count_cs = cs_max_pip_eqtl > 0.9,
  count_cs_99 = cs_99_max_pip_eqtl > 0.9
) %$%
  table(prediction, count_lead, useNA = "ifany") %>%
  as.matrix() %>%
  epitools::riskratio()

df.gbmi %>%
  dplyr::filter(n_r2 > 1) %$%
  table(prediction, (max_neff_r2 / min_neff_r2) > 2) %>%
  as.matrix() %>%
  epitools::riskratio()

df.gbmi %>%
  dplyr::filter(n_r2 > 1) %>%
  dplyr::group_by(prediction) %>%
  dplyr::summarize(median(max_neff_r2 / min_neff_r2))

########################################
trait <- "COPD_Bothsex_inv_var_meta_GBMI_052021"
region <- "chr1:161030340-162030340"


variant1 <- "chr1:161530340:A:G"
variant2 <- "chr1:161544752:A:C"
lead_pos1 <- 161530340
lead_pos2 <- 161544752

data <- munge_gbmi_sumstats(trait, region)
lead_idx1 <- which(data$variant == variant1)
lead_idx2 <- which(data$variant == variant2)
lead_rsid1 <- data$rsid[lead_idx1]
lead_rsid2 <- "rs396991"

plts <- plot_locuszoom_slalom(
  data,
  trait,
  region,
  additional_lead_pos = lead_pos2,
  window_size = 30000,
  fontsize = 8
)
p1 <- plts[[1]]
p2 <- plts[[2]]

p1[[1]] <- p1[[1]] +
  geom_text(
    aes(
      x = lead_pos2,
      y = -log10(data$p[lead_idx2]),
      label = "rs396991"
    ),
    data = tibble::tibble(),
    size = 2,
    vjust = -1
  ) +
  labs(tag = "d")
p1[[3]] <- p1[[3]] + guides(color = guide_legend(ncol = 3))
p2 <- p2 +
  geom_text(
    aes(
      x = data$r2[c(lead_idx1, lead_idx2)],
      y = data$chisq[c(lead_idx1, lead_idx2)],
      label = c(lead_rsid1, lead_rsid2)
    ),
    data = tibble::tibble(),
    size = 2,
    vjust = -1
  ) +
  labs(tag = "e") +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(hjust = 0.05))

##############################
my_pal <- function(range = c(1, 6)) {
  force(range)
  function(x) {
    scales::rescale(x, to = range, from = c(0, 1))
  }
}

df.ind <- rgsutil::read_gsfile(
  "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/COPD_Bothsex_inv_var_meta_GBMI_052021.chr1.161030340-162030340.full.txt.bgz"
) %>%
  dplyr::mutate(variant = stringr::str_c("chr", SNP)) %>%
  dplyr::filter(variant %in% (dplyr::filter(data, r2 > 0.6) %>% dplyr::pull(variant)))

df.N <- read.table(
  "./metadata/gbmi_sample_size_by_ancestry.txt",
  T,
  sep = "\t"
) %>%
  dplyr::filter(Endpoint == "COPD") %>%
  dplyr::mutate(
    cohort = ifelse(X..Total == "Total", "Meta", Biobank),
    ancestry = dplyr::case_when(
      cohort == "Meta" ~ "Meta",
      TRUE ~ stringr::str_to_upper(Ancestry)
    ),
    Ncase = Case,
    Ncontrol = Control,
    N = Ncase + Ncontrol,
    phi = Ncase / N,
    Neff = N * phi * (1 - phi)
  ) %>%
  dplyr::select(cohort, ancestry, N, Ncase, Ncontrol, Neff)

df.z <- dplyr::select(
  df.ind,
  variant,
  dplyr::ends_with("_beta"),
  dplyr::ends_with("_sebeta")
) %>%
  tidyr::pivot_longer(
    -variant,
    names_to = c("cohort", "ancestry", "variable"),
    names_pattern = "(.+)_(.+)_(.+)"
  ) %>%
  dplyr::filter(variant %in% c(variant1, variant2) &
    cohort != "all") %>%
  dplyr::mutate(
    variant = dplyr::case_when(
      variant == variant1 ~ "variant1",
      variant == variant2 ~ "variant2"
    ),
    cohort = ifelse(cohort == "all_inv_var", "Meta", cohort),
    ancestry = ifelse(cohort == "Meta", "Meta", stringr::str_to_upper(ancestry))
  ) %>%
  tidyr::pivot_wider(
    names_from = c("variant", "variable"),
    values_from = "value"
  ) %>%
  dplyr::left_join(df.N) %>%
  dplyr::mutate(
    missingness = dplyr::case_when(
      is.na(variant1_beta) ~ sprintf("Missing %s", lead_rsid1),
      is.na(variant2_beta) ~ sprintf("Missing %s", lead_rsid2),
      TRUE ~ "Both exist"
    ),
    variant1_z = ifelse(is.na(variant1_beta), 0, variant1_beta / variant1_sebeta),
    variant2_z = ifelse(is.na(variant2_beta), 0, variant2_beta / variant2_sebeta),
    variant1_beta = ifelse(is.na(variant1_beta), 8, variant1_beta),
    variant2_beta = ifelse(is.na(variant2_beta), 8, variant2_beta),
    label = dplyr::case_when(
      missingness != "Both exist" ~ cohort,
      cohort == "Meta" ~ cohort,
      # abs(variant1_z) > 3 ~ cohort,
      # abs(variant2_z) > 3 ~ cohort,
      TRUE ~ ""
    ),
    ancestry = factor(ancestry, levels = c(sort(
      setdiff(unique(ancestry), "Meta")
    ), "Meta"))
  )

dplyr::filter(df.z, !is.na(variant2_sebeta) & cohort != "Meta") %>%
  dplyr::left_join(df.N) %>%
  dplyr::summarize(
    Ncase = sum(Ncase),
    Ncontrol = sum(Ncontrol),
    N = Ncase + Ncontrol,
    phi = Ncase / N,
    Neff = N * phi * (1 - phi)
  )

p3 <-
  ggplot(df.z, aes(variant1_z, variant2_z, color = ancestry)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(aes(size = Neff, shape = missingness)) +
  ggrepel::geom_text_repel(
    aes(
      label = label,
      point.size = Neff
    ),
    size = 2,
    segment.size = unit(0.25, "mm"),
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
    box.padding = unit(0.34, "lines"),
    seed = 1000
  ) +
  locusviz::get_default_theme(legend.position = "right") +
  theme(legend.box = "horizontal") +
  coord_cartesian(xlim = c(-2.5, 7), ylim = c(-2.5, 7)) +
  scale_shape_manual(values = c(21, 16)) +
  scale_color_manual(
    values = gnomad_pop_colors,
    guide = guide_legend(override.aes = list(label = ""), order = 1)
  ) +
  continuous_scale(
    aesthetics = c("size", "point.size"),
    scale_name = "size",
    breaks = c(0, 1000, 5000, 10000, 50000, 100000),
    labels = scales::comma_format(),
    palette = my_pal(c(1, 3)),
    guide = guide_legend(
      title = "Effective N",
      override.aes = list(label = "")
    )
  ) +
  labs(
    x = sprintf("Z-score (%s)", lead_rsid1),
    y = sprintf("Z-score (%s)", lead_rsid2),
    color = "Ancestry",
    shape = "Missingness",
    tag = "f"
  )


layout <- "
ABC
DDD
EEE
FFF
GGG
HI#
"

plt <- purrr::reduce(c(
  plts_validation,
  p1,
  list(
    p2,
    p3 + theme(
      legend.position = "none",
      plot.tag.position = c(0.05, 1)
    )
  )
), `+`) + patchwork::plot_layout(
  heights = c(1.4, 1, 0.35, 0.35, 0.05, 1.4),
  design = layout
)
cowplot::save_plot(
  "figures/Fig6/Fig6_gbmi_validation_base.pdf",
  plt,
  base_height = 7.4,
  base_width = 7.2,
  device = cairo_pdf
)

cowplot::save_plot(
  "figures/Fig6/Fig6_beta_plot_legend.pdf",
  p3,
  base_height = 2.4,
  base_width = 4.8,
  device = cairo_pdf
)