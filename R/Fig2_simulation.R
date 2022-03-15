library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

df <- dplyr::bind_rows(
  rgsutil::read_gsfile(
    "gs://meta-finemapping-simulation/new_simulations/sim.ALL.metrics.tsv.bgz"
  ),
  rgsutil::read_gsfile(
    "gs://meta-finemapping-simulation/new_simulations/sim_grch38.ALL.metrics.tsv.bgz"
  ) %>%
    dplyr::filter(pheno_group %in% c(1, 3, 5)) %>%
    dplyr::mutate(config_group = -config_group)
)

parse_array <- function(x) {
  as.numeric(t(stringr::str_split_fixed(
    stringr::str_remove_all(x, "[\\[\\]]"), ",", 5
  )))
}

pip_bin_breaks <- c(0, 0.01, 0.1, 0.5, 0.9, 1.0)
q_bin_breaks <- c(0, 1e-4, 5e-4, 1e-3, 5e-3, 0.01, 1)

get_metrics <- function(df) {
  df %>%
    tidyr::nest(-pheno_group, -config_group, -Ngamma_region) %>%
    dplyr::mutate(data = purrr::map(data, function(x) {
      data.frame(
        prob_bin = pip_bin_breaks[1:5],
        q_bin = q_bin_breaks[2:6],
        calibration = parse_array(x$calibration),
        TP = parse_array(x$prob_true_bins),
        Nbin = parse_array(x$prob_bins),
        mean_prob = parse_array(x$mean_prob),
        mean_prob_sd = parse_array(x$mean_prob_sd),
        q_true_bins = parse_array(x$q_true_bins)
      ) %>%
        dplyr::mutate(
          Ngamma = sum(TP),
          TP.q.cumsum = cumsum(tidyr::replace_na(q_true_bins, 0))
        )
    })) %>%
    tidyr::unnest(data) %>%
    dplyr::mutate(
      pheno_group = as.factor(pheno_group),
      prob_bin = factor(prob_bin),
      calibration.lower = na_binconf(TP, Nbin, "Lower"),
      calibration.upper = na_binconf(TP, Nbin, "Upper"),
      fdr = 1 - calibration,
      fdr.lower = 1 - calibration.upper,
      fdr.upper = 1 - calibration.lower,
      n_mean_prob = 1 - mean_prob,
      q_bin = factor(q_bin, labels = 100 * unique(q_bin)),
      recall = TP.q.cumsum / Ngamma_region,
      recall.lower = na_binconf(TP.q.cumsum, Ngamma_region, "Lower"),
      recall.upper = na_binconf(TP.q.cumsum, Ngamma_region, "Upper")
    )
}

metrics <- get_metrics(df) %>%
  dplyr::mutate(
    array = factor(
      dplyr::case_when(
        config_group %in% c(-1, seq(1, 3), seq(10, 12), 1822, 3337) ~ "Omni25",
        config_group %in% seq(4, 6) ~ "MEGA",
        config_group %in% seq(7, 9) ~ "GSA",
        config_group %in% c(1317, 2327, 2832, 3842) ~ "Mixed",
        TRUE ~ as.character(config_group)
      ),
      levels = c("Omni25", "MEGA", "GSA", "Mixed")
    ),
    imputation = factor(
      dplyr::case_when(
        config_group %in% c(-1) ~ "TopMed",
        config_group %in% c(1, 4, 7, seq(10, 12), 1317, 2832) ~ "1000GP3",
        config_group %in% c(2, 5, 8) ~ "HRC",
        config_group %in% c(3, 6, 9) ~ "TopMed-liftover",
        config_group %in% c(1822, 2327, 3337, 3842) ~ "Mixed",
        TRUE ~ as.character(config_group)
      ),
      levels = c("1000GP3", "HRC", "TopMed", "TopMed-liftover", "Mixed")
    ),
    ancestry = factor(
      dplyr::case_when(
        config_group %in% c(10) ~ "AFR + EUR",
        config_group %in% c(11) ~ "EAS + EUR",
        config_group %in% c(12, 2832, 3337, 3842) ~ "AFR + EAS + EUR",
        TRUE ~ "EUR"
      ),
      levels = c(
        "EUR",
        "EAS + EUR",
        "AFR + EUR",
        "AFR + EAS + EUR"
      )
    ),
    mixture = factor(
      dplyr::case_when(
        config_group == 2327 ~ "Mixed array & panel (EUR)",
        config_group == 3842 ~ "Mixed array & panel (multi-ancestry)",
        TRUE ~ NA_character_
      ),
      levels = c(
        "Mixed array & panel (EUR)",
        "Mixed array & panel (multi-ancestry)"
      )
    )
  )

##################################################

dplyr::filter(metrics, !is.na(mixture) & prob_bin == 0.9) %>%
  dplyr::mutate(expected = 1 - mean_prob) %>%
  View()
dplyr::filter(
  metrics,
  imputation == "1000GP3" &
    ancestry == "EUR" & prob_bin == 0.9
) %>%
  dplyr::mutate(expected = 1 - mean_prob) %>%
  View()

dplyr::filter(
  metrics,
  array == "Omni25" &
    ancestry == "EUR" & prob_bin == 0.9
) %>%
  dplyr::mutate(expected = 1 - mean_prob) %>%
  View()

##################################################
plot_fdr_and_recall <- function(df,
                                x_var,
                                xlab,
                                cols,
                                no_y = FALSE,
                                no_shape_legend = TRUE,
                                prob_bin = 0.9,
                                q_bin = 0.01) {
  df1 <- dplyr::filter(df, prob_bin == .env$prob_bin)
  print(levels(df1[[x_var]]))
  df2 <- dplyr::filter(df, q_bin == .env$q_bin)
  levels(df2[[x_var]]) <- stringr::str_replace(levels(df2[[x_var]]), " \\(", "\n\\(")
  position <- position_dodge(width = 0.75)

  p1 <-
    ggplot(
      df1,
      aes_string(
        x_var,
        "fdr",
        color = x_var,
        group = "pheno_group",
        shape = "pheno_group"
      )
    ) +
    geom_point(
      aes_string(x_var, "n_mean_prob"),
      color = "grey50",
      size = 3,
      shape = "-",
      position = position
    ) +
    geom_errorbar(aes(ymin = fdr.lower, ymax = fdr.upper),
      width = 0,
      position = position
    ) +
    geom_point(size = 1, position = position) +
    scale_color_manual(values = c(cols, NA), drop = FALSE) +
    locusviz::or_missing(no_shape_legend, scale_shape(guide = FALSE)) +
    locusviz::or_missing(
      !no_shape_legend,
      guides(
        color = guide_legend(order = 1),
        shape = guide_legend(
          order = 2,
          direction = "horizontal",
          title.vjust = 1,
          label.vjust = 1
        )
      )
    ) +
    labs(
      x = xlab,
      y = "False discovery rate",
      color = xlab,
      shape = expression("Effect size " ~ italic(r)[g])
    ) +
    coord_cartesian(ylim = c(0, 0.6)) +
    locusviz::get_default_theme(hide.ylab = no_y) +
    theme(
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.margin = margin(2, 4, -2, 4, unit = "pt")
    )

  p2 <-
    ggplot(
      df2,
      aes_string(
        x_var,
        "recall",
        color = x_var,
        group = "pheno_group",
        shape = "pheno_group"
      )
    ) +
    geom_errorbar(
      aes(ymin = recall.lower, ymax = recall.upper),
      width = 0,
      position = position
    ) +
    geom_point(size = 1, position = position) +
    scale_color_manual(values = cols) +
    scale_shape(guide = FALSE) +
    labs(
      x = xlab,
      y = "Recall",
      color = xlab
    ) +
    coord_cartesian(ylim = c(0, 0.6)) +
    locusviz::get_default_theme(hide.ylab = no_y) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      )
    )
  return(list(p1, p2))
}


plt <-
  list(
    dplyr::filter(
      metrics,
      imputation == "1000GP3" & ancestry == "EUR"
    ) %>%
      dplyr::mutate(pheno_group = factor(
        dplyr::case_when(
          pheno_group == 1 ~ 1,
          pheno_group == 3 ~ 0.9,
          pheno_group == 5 ~ 0.5
        ),
        levels = c(1, 0.9, 0.5)
      )) %>%
      plot_fdr_and_recall(
        x_var = "array",
        xlab = "Genotyping array",
        no_shape_legend = FALSE,
        cols = c(
          BuenColors::jdb_palette("calma_morado")[c(4, 5, 7)],
          BuenColors::jdb_palette("brewer_blue")[9]
        )
      ),
    dplyr::filter(
      metrics,
      array == "Omni25" & ancestry == "EUR"
    ) %>%
      dplyr::mutate(imputation = factor(
        dplyr::case_when(
          imputation == "1000GP3" ~ "1000GP3 (matching ref)",
          TRUE ~ as.character(imputation)
        ),
        levels = c(
          "1000GP3 (matching ref)",
          "HRC",
          "TopMed",
          "TopMed-liftover",
          "Mixed"
        )
      )) %>%
      plot_fdr_and_recall(
        x_var = "imputation",
        xlab = "Imputation panel",
        cols = c(
          BuenColors::jdb_palette("brewer_marine")[7],
          BuenColors::jdb_palette("GrandBudapest")[2],
          BuenColors::jdb_palette("Rushmore")[3],
          BuenColors::jdb_palette("brewer_marine")[5],
          BuenColors::jdb_palette("brewer_blue")[9]
        ),
        no_y = TRUE
      ),
    dplyr::filter(
      metrics,
      array == "Omni25" & imputation == "1000GP3"
    ) %>%
      plot_fdr_and_recall(
        x_var = "ancestry",
        xlab = "Genetic ancestry",
        cols = c(
          BuenColors::jdb_palette("brewer_blue")[8],
          BuenColors::jdb_palette("brewer_green")[7],
          BuenColors::jdb_palette("brewer_fire")[c(6, 8)]
        ),
        no_y = TRUE
      ),
    dplyr::filter(metrics, !is.na(mixture)) %>%
      plot_fdr_and_recall(
        x_var = "mixture",
        xlab = "Heterogenous settings",
        cols = c(
          BuenColors::jdb_palette("brewer_blue")[9],
          BuenColors::jdb_palette("brewer_orange")[9]
        ),
        no_y = TRUE
      )
  ) %>%
  purrr::flatten() %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 2, byrow = FALSE) +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "figures/Fig2/Fig2_simulation.pdf",
  plt,
  base_height = 4.0,
  base_width = 7.2,
  device = cairo_pdf
)
