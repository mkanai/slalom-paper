library(ggplot2)
library(dplyr)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

levels <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/imputed_info/level.txt") %>%
  dplyr::mutate(label = stringr::str_replace(label, "TopMed", "TOPMed-liftover"))

plot_upset <- function(maf_str, rsq_str, maf_prefix = "") {
  df.count <- rgsutil::read_gsfile(
    sprintf(
      "gs://meta-finemapping-simulation/imputed_info/upset.%smaf%s.rsq%s.tsv.bgz",
      maf_prefix,
      maf_str,
      rsq_str
    )
  )

  print(dplyr::summarize(
    df.count,
    numer = n[which(collection == "abcdefghijklmnopqrstu")],
    denom = sum(n),
    frac = numer / denom
  ))

  mat <-
    dplyr::arrange(df.count, -n) %>%
    head(30) %>%
    dplyr::mutate(level = stringr::str_split(collection, "")) %>%
    tidyr::unnest(level) %>%
    dplyr::left_join(levels) %>%
    tidyr::separate(
      label,
      sep = "_",
      into = c("ancestry", "array", "imputation"),
      remove = FALSE
    ) %>%
    dplyr::mutate(
      ancestry = factor(ancestry, levels = c("AFR", "EAS", "EUR")),
      array = factor(array, levels = c("GSA", "MEGA", "Omni25")),
      imputation = factor(imputation, levels = c("1000GP3", "HRC", "TOPMed-liftover")),
    ) %>%
    dplyr::group_by(collection) %>%
    dplyr::mutate(n_pop = length(unique(ancestry))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      collection = factor(collection, levels = unique(collection[order(-n_pop, -n)])),
      label = factor(label, levels = unique(label[order(ancestry, imputation, array)]), ordered = TRUE)
    )

  segment <-
    dplyr::group_by(mat, collection) %>%
    dplyr::summarize(y = min(label), yend = max(label))

  bar <- dplyr::group_by(mat, collection) %>%
    dplyr::summarize(
      n = n[1],
      ancestry = length(unique(ancestry)),
      array = length(unique(array)),
      imputation = length(unique(imputation))
    )

  p_mat <-
    ggplot() +
    geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = fill
      ),
      data = tibble::tibble(
        xmin = -Inf,
        xmax = Inf,
        ymin = -0.5 + seq(21),
        ymax = 0.5 + seq(21)
      ),
      fill = c(rep("#F2F2F2", 9), rep("white", 6), rep("#F2F2F2", 6)),
    ) +
    geom_hline(
      yintercept = 0.5 + seq(3, 18, by = 3),
      size = 0.1,
      color = "grey50"
    ) +
    geom_vline(
      xintercept = 0.5 + cumsum(purrr::map_dbl(3:2, function(x) {
        sum(bar$ancestry == x)
      })),
      linetype = "dashed",
      color = "grey50"
    ) +
    scale_x_discrete() +
    scale_y_discrete(limits = rev(levels(mat$label))) +
    geom_segment(
      aes(
        x = collection,
        xend = collection,
        y = y,
        yend = yend
      ),
      size = 0.1,
      color = "grey50",
      data = segment
    ) +
    geom_point(aes(collection, label, color = imputation, shape = array), data = mat) +
    locusviz::get_default_theme(hide.xlab = TRUE, hide.ytitle = TRUE) +
    theme(
      legend.position = "right",
      legend.justification = c(0, 1),
      axis.ticks.x = element_blank(),
      legend.margin = margin()
    ) +
    labs(color = "Imputation\npanel", shape = "Genotyping\narray") +
    scale_color_manual(
      values = c(
        BuenColors::jdb_palette("brewer_marine")[7],
        BuenColors::jdb_palette("GrandBudapest")[2],
        BuenColors::jdb_palette("Rushmore")[3]
      )
    ) +
    guides(color = guide_legend(order = 1))

  p_bar <-
    ggplot() +
    geom_vline(
      xintercept = 0.5 + cumsum(purrr::map_dbl(3:2, function(x) {
        sum(bar$ancestry == x)
      })),
      linetype = "dashed",
      color = "grey50"
    ) +
    scale_x_discrete() +
    geom_point(aes(
      collection,
      n,
      color = factor(imputation),
      shape = factor(array)
    ), data = bar) +
    # geom_text(
    #   aes(collection, n, label = scales::label_number_si()(n)),
    #   size = 2,
    #   vjust = -1,
    #   data = bar
    # ) +
    ggrepel::geom_text_repel(
      aes(collection, n, label = scales::label_number_si()(n)),
      size = 2,
      point.size = 1,
      segment.size = 0,
      nudge_y = 0.2,
      box.padding = unit(0.02, "line"),
      direction = "y",
      data = bar
    ) +
    scale_y_log10(
      label = scales::comma_format(),
      expand = expansion(0.2, 0)
    ) +
    scale_shape_manual(values = c(21, 24, 22)) +
    scale_color_manual(
      values = c(
        BuenColors::jdb_palette("brewer_purple")[7],
        BuenColors::jdb_palette("Darjeeling")[4],
        BuenColors::jdb_palette("Darjeeling")[1]
      ),
      guide = guide_legend(override.aes = aes(shape = 21))
    ) +
    locusviz::get_default_theme(hide.xlab = TRUE) +
    theme(
      legend.position = "right",
      legend.justification = c(0, 1),
      legend.margin = margin()
    ) +
    labs(y = "# variants", shape = "# genotyping\n   arrays", color = "# imputation\n   panels")

  return(list(
    p_bar,
    p_mat
  ))
}

plts <- purrr::pmap(tibble::tibble(
  maf_str = c("0001", "005", "0001", "0001"),
  rsq_str = "06",
  maf_prefix = c("gnomad_min_", "gnomad_min_", "gnomad_max_min_", ""),
  i = seq(4)
), function(maf_str, rsq_str, maf_prefix, i) {
  plts <- plot_upset(maf_str, rsq_str, maf_prefix)
  title <- sprintf(
    "%s MAF > %s%s & Rsq > 0.6",
    ifelse(
      stringr::str_starts(maf_prefix, "gnomad"),
      "gnomAD",
      "per-combination"
    ),
    stringr::str_replace(maf_str, "^0", "0."),
    dplyr::case_when(
      maf_prefix == "gnomad_min_" ~ " in every ancestry",
      maf_prefix == "gnomad_max_min_" ~ " in >1 ancesty but not every",
      TRUE ~ ""
    ),
    stringr::str_replace(rsq_str, "^0", "0.")
  )
  plts[[1]] <- plts[[1]] + labs(
    tag = letters[2 - i %% 2],
    title = title
  ) +
    theme(plot.title = element_text(margin = margin(t = -7), hjust = -0.013 * stringr::str_length(title)))
  return(plts)
})

p1 <-
  magrittr::extract(plts, c(1, 2)) %>%
  purrr::flatten() %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1, heights = rep(c(1, 2), 2))

p2 <-
  magrittr::extract(plts, c(3, 4)) %>%
  purrr::flatten() %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 1, heights = rep(c(1, 2), 2))

cowplot::save_plot(
  "figures/SFig3_imputation_upset.pdf",
  p1,
  base_height = 8,
  base_width = 7.2,
  device = cairo_pdf
)
cowplot::save_plot(
  "figures/SFig4_imputation_upset.pdf",
  p2,
  base_height = 8,
  base_width = 7.2,
  device = cairo_pdf
)
