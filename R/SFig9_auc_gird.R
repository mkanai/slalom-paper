library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

future::plan("multisession")
df.auc.grid <-
  furrr::future_pmap_dfr(tidyr::crossing(
    r2_threshold = seq(0, 0.9, by = 0.1),
    nlog10p_dentist_lower_bound = c(seq(0, 10, 0.5))
  ), function(r2_threshold,
              nlog10p_dentist_lower_bound) {
    df <- dplyr::filter(df.train, r2_bin >= r2_threshold &
      max_prob_bin >= 0.9) %>%
      munge_metrics(nlog10p_dentist_lower_bound = nlog10p_dentist_lower_bound)

    truth <- c("abf_success", "abf_success_cs", "abf_success_cs_99")
    roc_auc <- purrr::map_dbl(truth, function(x) {
      yardstick::roc_auc(df, !!as.symbol(x), "count")$.estimate
    })
    tibble::tibble(
      r2_threshold = r2_threshold,
      nlog10p_dentist_lower_bound = nlog10p_dentist_lower_bound,
      truth = truth,
      roc_auc = roc_auc
    )
  })


plt <- purrr::map(c("abf_success", "abf_success_cs", "abf_success_cs_99"), function(truth) {
  label <- dplyr::case_when(
    truth == "abf_success" ~ "Lead PIP variant",
    truth == "abf_success_cs" ~ "in 95% CS",
    truth == "abf_success_cs_99" ~ "in 99% CS",
  )
  dplyr::filter(df.auc.grid, truth == .env$truth) %>%
    ggplot(aes(r2_threshold, nlog10p_dentist_lower_bound)) +
    geom_tile(aes(fill = roc_auc)) +
    geom_rect(
      aes(
        xmin = 0.55,
        xmax = 0.65,
        ymin = 3.75,
        ymax = 4.25
      ),
      color = "grey50",
      fill = NA
    ) +
    scale_fill_viridis_c() +
    labs(
      x = expression(paste(italic(r)^2, " threshold")),
      y = expression(paste(-log[10], " ", italic(P)[`DENTIST-S`], " threshold")),
      fill = "AUROC",
      title = label
    ) +
    locusviz::get_default_theme(
      legend.position = "bottom",
      legend.justification = "bottom",
      hide.ylab = truth != "abf_success"
    ) +
    theme(
      legend.direction = "horizontal",
      legend.key.width = unit(0.8, "cm"),
      plot.title = element_text(hjust = 2e-2),
    ) +
    scale_x_continuous(expand = expansion()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
}) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/SFig9_grid_auc.pdf",
  plt,
  base_height = 3.6,
  base_width = 7.2,
  device = cairo_pdf
)