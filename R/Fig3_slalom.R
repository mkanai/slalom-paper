library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

########################################################

generate_example_data <- function(lead_chisq = 50, seed = 123456) {
  set.seed(seed)
  r2 <- c(
    runif(3, 0, 0.05),
    runif(5, 0, 0.3),
    0.95,
    0.99,
    1,
    0.96,
    0.93,
    0.9,
    runif(5, 0, 0.3)
  )
  chisq_exp <- lead_chisq * r2
  chisq_obs <- pmax(chisq_exp - rnorm(length(chisq_exp)), 1)
  chisq_obs[1:3] <- lead_chisq * runif(3, 0.65, 0.85)
  chisq_obs[11] <- lead_chisq
  chisq_obs[12] <- chisq_obs[12] * 0.65
  shape <- rep(16, length(chisq_obs))
  shape[11] <- 18
  shape[12] <- 17

  data <- tibble::tibble(
    chisq = chisq_obs,
    x = seq(length(chisq)),
    nlog10p = -log10(pchisq(chisq_obs, 1, lower.tail = F)),
    r2 = r2,
    shape = factor(shape)
  )
  return(data)
}

plot_example_locuszoom <- function(data) {
  plt <-
    dplyr::mutate(data, r2 = ifelse(x == 11, NA, r2)) %>%
    ggplot(aes(x, nlog10p)) +
    geom_hline(
      yintercept = -log10(5e-8),
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_point(aes(
      color = r2,
      shape = shape,
      size = shape
    )) +
    locusviz::get_default_theme(hide.xtext = TRUE) +
    theme(axis.ticks.x = element_blank()) +
    labs(
      x = "Locus",
      y = expression(paste(-log[10], "(", italic(P), ")")),
      color = expression(italic(r)^2)
    ) +
    scale_color_stepsn(
      colors = c("navy", "lightskyblue", "green", "orange", "red"),
      breaks = seq(0.2, 0.8, by = 0.2),
      limits = c(0, 1),
      show.limits = TRUE,
      na.value = "purple3"
    ) +
    scale_shape_manual(
      values = c(
        `16` = 16,
        `17` = 18,
        `18` = 18
      ),
      guide = FALSE
    ) +
    scale_size_manual(
      values = c(
        `16` = 1.5,
        `17` = 3,
        `18` = 3
      ),
      guide = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, 12.5),
      breaks = scales::pretty_breaks(),
      expand = expansion(0, 0)
    )

  return(plt)
}

plot_dentist_distribution <- function(lead_chisq,
                                      nlog10p_dentist_lower_bound = 4,
                                      legend = FALSE,
                                      n_grid = 1000) {
  data <- expand.grid(
    r2 = seq(0, 1, length.out = n_grid),
    chisq = seq(0, lead_chisq, length.out = n_grid)
  ) %>%
    dplyr::mutate(
      t = (sqrt(chisq) - sqrt(r2 * lead_chisq))^2 / (1 - r2),
      nlog10p_dentist = pchisq(t, 1, lower.tail = F, log.p = TRUE) / -log(10)
    )

  ggplot(data, aes(r2, chisq)) +
    geom_raster(aes(fill = nlog10p_dentist), interpolate = TRUE) +
    locusviz::or_missing(
      lead_chisq > qchisq(5e-8, 1, lower.tail = F),
      geom_hline(
        yintercept = qchisq(5e-8, 1, lower.tail = F),
        linetype = "dashed",
        color = "white"
      )
    ) +
    geom_abline(
      intercept = 0,
      slope = lead_chisq,
      linetype = "dashed",
      color = "white"
    ) +
    scale_fill_gradientn(
      colors = BuenColors::jdb_palette("solar_extra"),
      limits = c(0, nlog10p_dentist_lower_bound),
      oob = scales::squish
    ) +
    scale_x_continuous(
      expand = expansion(0, 0),
      labels = scales::number_format(drop0trailing = TRUE)
    ) +
    scale_y_continuous(
      expand = expansion(0, 0),
      labels = scales::comma_format(drop0trailing = TRUE)
    ) +
    locusviz::get_default_theme(fontsize = 8) +
    theme(plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm")) +
    coord_fixed(ratio = 1 / lead_chisq, clip = "off") +
    locusviz::or_missing(
      !legend,
      theme(legend.position = "none")
    ) +
    locusviz::or_missing(
      legend,
      theme(
        legend.position = "bottom",
        legend.key.width = unit(0.4, "cm")
      )
    ) +
    labs(
      x = expression(paste(italic(r)^2, "to the lead variant")),
      y = expression(italic(chi)^2),
      fill = expression(paste(-log[10], " ", italic(P)[`DENTIST-S`]))
    )
}

data <- generate_example_data()
p1 <- plot_example_locuszoom(data)
p2 <- plot_dentist_distribution(lead_chisq = 50, legend = TRUE) +
  ggnewscale::new_scale_fill() +
  geom_point(
    aes(
      r2,
      chisq,
      shape = shape,
      size = shape,
      fill = shape
    ),
    color = "grey20",
    data = data
  ) +
  scale_shape_manual(
    values = c(`16` = 21, `17` = 23, `18` = 23),
    guide = FALSE
  ) +
  scale_size_manual(values = c(`16` = 1.5, `17` = 3, `18` = 3), guide = FALSE) +
  scale_fill_manual(
    values = c(`16` = "grey70", `17` = "white", `18` = "purple3"),
    guide = FALSE
  )

########################################################
plt_grid <-
  purrr::map(c(10, 50, 100, 500), function(lead_chisq) {
    plot_dentist_distribution(lead_chisq = lead_chisq, legend = (lead_chisq == 500))
  }) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 2) +
  patchwork::plot_annotation(tag_levels = "a")
plt_grid


########################################################
df.frac <-
  dplyr::bind_rows(
    rgsutil::read_gsfile(
      "gs://meta-finemapping-simulation/new_simulations/sim.ALL.dentist.frac.tsv.bgz"
    ),
    rgsutil::read_gsfile(
      "gs://meta-finemapping-simulation/new_simulations/sim_grch38.ALL.dentist.frac.tsv.bgz"
    )
  ) %>%
  tidyr::drop_na(r2_bin, nlog10p_dentist_bin) %>%
  dplyr::mutate(
    r2_bin = r2_bin * 0.1,
    max_prob_bin = max_prob_bin * 0.1,
    nlog10p_dentist_bin = nlog10p_dentist_bin * 0.5
  )


munge_metrics <- function(df, nlog10p_dentist_lower_bound = 4) {
  dplyr::group_by(df, pheno, config, region) %>%
    dplyr::summarize(
      abf_success = factor(any(abf_success), levels = c("FALSE", "TRUE")),
      abf_success_cs = factor(any(abf_success_cs), levels = c("FALSE", "TRUE")),
      abf_success_cs_99 = factor(any(abf_success_cs_99), levels = c("FALSE", "TRUE")),
      max_prob_bin = max_prob_bin[1],
      max_nlog10p_dentist = max(max_nlog10p_dentist, na.rm = T),
      count = as.numeric(sum((nlog10p_dentist_bin >= nlog10p_dentist_lower_bound) * n)),
      count_total = as.numeric(sum(n)),
      fraction = count / count_total,
    ) %>%
    dplyr::ungroup()
}


df.train <- df.frac

df.train.best <-
  dplyr::filter(df.train, r2_bin >= r2_threshold &
    max_prob_bin >= 0.9) %>%
  munge_metrics()

df.train.prediction <- dplyr::mutate(
  dplyr::filter(df.train, r2_bin >= r2_threshold) %>%
    munge_metrics(),
  prediction = factor(
    count >= count_threshold,
    levels = c(TRUE, FALSE),
    labels = c("Suspicious loci", "Non-suspicious loci")
  )
)

df.prob001 <-
  dplyr::bind_rows(
    rgsutil::read_gsfile(
      "gs://meta-finemapping-simulation/new_simulations/sim.ALL.prob001.tsv.bgz"
    ),
    rgsutil::read_gsfile(
      "gs://meta-finemapping-simulation/new_simulations/sim_grch38.ALL.prob001.tsv.bgz"
    ),
  ) %>%
  dplyr::left_join(dplyr::select(df.train.prediction, pheno, config, region, prediction)) %>%
  dplyr::mutate(prob_bin = cut(prob, pip_bin_breaks))

df.prob001.threshold <-
  purrr::map_dfr(seq(0.01, 1, 0.01), function(threshold) {
    dplyr::filter(df.prob001, prob > threshold & !is.na(prediction)) %>%
      dplyr::group_by(prediction) %>%
      dplyr::summarize(
        threshold = threshold,
        mean_prob = mean(prob, na.rm = T),
        calibration = mean(gamma),
        calibration.lower = na_binconf(sum(gamma), n(), "Lower"),
        calibration.upper = na_binconf(sum(gamma), n(), "Upper"),
        fdr = 1 - calibration,
        fdr.lower = 1 - calibration.upper,
        fdr.upper = 1 - calibration.lower,
        count = n()
      )
  })

plot_curves <- function(df, predictor, threshold) {
  truth <- c("abf_success", "abf_success_cs", "abf_success_cs_99")
  roc_auc <- purrr::map_dbl(truth, function(x) {
    yardstick::roc_auc(df, !!as.symbol(x), !!as.symbol(predictor))$.estimate
  })
  sens <- purrr::map_dbl(truth, function(x) {
    yardstick::sensitivity(df, !!as.symbol(x), factor(!!as.symbol(predictor) < threshold))$.estimate
  })
  spec <- purrr::map_dbl(truth, function(x) {
    yardstick::specificity(df, !!as.symbol(x), factor(!!as.symbol(predictor) < threshold))$.estimate
  })
  prec <- purrr::map_dbl(truth, function(x) {
    yardstick::precision(df, !!as.symbol(x), factor(!!as.symbol(predictor) < threshold))$.estimate
  })

  df.label <- tibble::tibble(
    x = 1 - spec,
    y = sens,
    truth = truth,
    label = ifelse(truth == "abf_success_cs_99", "# outlier variants > 0", "")
  )

  p1 <- dplyr::bind_rows(
    yardstick::roc_curve(df, abf_success, !!as.symbol(predictor)) %>%
      dplyr::mutate(
        truth = sprintf("Lead PIP variant\n(AUROC = %.2f)", roc_auc[1]),
        level = 1
      ),
    yardstick::roc_curve(df, abf_success_cs, !!as.symbol(predictor)) %>%
      dplyr::mutate(
        truth = sprintf("in 95%% CS\n(AUROC = %.2f)", roc_auc[2]),
        level = 2
      ),
    yardstick::roc_curve(df, abf_success_cs_99, !!as.symbol(predictor)) %>%
      dplyr::mutate(
        truth = sprintf("in 99%% CS\n(AUROC = %.2f)", roc_auc[3]),
        level = 3
      ),
  ) %>%
    dplyr::mutate(truth = factor(truth, levels = unique(truth[order(level)]))) %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_path(aes(color = truth)) +
    geom_point(aes(x, y),
      color = "black",
      size = 1,
      data = df.label
    ) +
    geom_text(
      aes(x, y, label = label),
      data = df.label,
      vjust = -1,
      hjust = 0.85,
      size = 2
    ) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      expand = expansion(0, 0),
      labels = scales::number_format(drop0trailing = TRUE)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(0, 0),
      labels = scales::number_format(drop0trailing = TRUE)
    ) +
    locusviz::get_default_theme() +
    theme(
      legend.position = c(1, 0.01),
      legend.justification = c(1, 0),
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_blank()
    ) +
    scale_color_manual(
      values = c(
        BuenColors::jdb_palette("Darjeeling")[4],
        BuenColors::jdb_palette("brewer_purple")[6],
        BuenColors::jdb_palette("Darjeeling")[2]
      )
    )
  return(p1)
}

plot_new_calibration <- function(data) {
  ggplot() +
    geom_rect(
      aes(
        xmin = 0,
        xmax = 0.1,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = "grey50",
      alpha = 0.5,
      data = tibble::tibble(
        xmin = 0,
        xmax = 0.1,
        ymin = -Inf,
        ymax = Inf
      )
    ) +
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
        mean_prob - calibration,
        ymin = mean_prob - calibration.lower,
        ymax = mean_prob - calibration.upper,
        fill = prediction
      ),
      alpha = 0.5,
      data = data
    ) +
    geom_path(aes(threshold, mean_prob - calibration, color = prediction),
      data = data
    ) +
    locusviz::get_default_theme() +
    scale_color_manual(values = BuenColors::jdb_palette("brewer_celsius")[c(8, 2)]) +
    scale_fill_manual(
      values = BuenColors::jdb_palette("brewer_celsius")[c(8, 2)],
      guide = "none"
    ) +
    labs(x = "PIP threshold", y = "Mean PIP - % true causal", color = "Prediction") +
    theme(
      legend.position = c(1, 1),
      legend.justification = c(1, 1)
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.1, 0.25, 0.5, 0.75, 1),
      expand = expansion(0, 0),
      labels = scales::number_format(drop0trailing = TRUE)
    ) +
    scale_y_continuous(
      limits = c(0, 0.35),
      labels = scales::number_format(drop0trailing = TRUE)
    )
}

plot_variant_fraction <- function(data) {
  dplyr::group_by(data, threshold) %>%
    dplyr::summarize(
      frac = count / sum(count),
      prediction = prediction
    ) %>%
    ggplot(aes(threshold, frac)) +
    geom_col(aes(fill = prediction),
      width = 0.01,
      position = "stack"
    ) +
    geom_rect(
      aes(
        x = NULL,
        y = NULL,
        xmin = 0,
        xmax = 0.1,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = "grey50",
      alpha = 0.5,
      data = tibble::tibble(
        xmin = 0,
        xmax = 0.1,
        ymin = -Inf,
        ymax = Inf
      )
    ) +
    scale_fill_manual(
      values = BuenColors::jdb_palette("brewer_celsius")[c(8, 2)],
      guide = "none"
    ) +
    locusviz::get_default_theme() +
    labs(x = "PIP threshold", y = "% variants") +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.1, 0.25, 0.5, 0.75, 1),
      expand = expansion(0, 0),
      labels = scales::number_format(drop0trailing = TRUE)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(0, 0),
      labels = scales::percent_format(drop0trailing = TRUE)
    )
}

p3 <- plot_curves(df.train.best, "count", count_threshold)
p4 <- plot_new_calibration(df.prob001.threshold)
p5 <- plot_variant_fraction(df.prob001.threshold)

plt <- list(p1, p2, patchwork::plot_spacer(), p3, p4, p5) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 2) +
  patchwork::plot_annotation(tag_levels = "a")
plt

cowplot::save_plot(
  "figures/Fig3/Fig3_slalom_overview.pdf",
  plt,
  base_height = 5.2,
  base_width = 7.2,
  device = cairo_pdf
)