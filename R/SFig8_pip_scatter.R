library(ggplot2)
library(dplyr)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

df.delta <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi.bbj_fg_ukbb.delta.tsv.bgz")
df.delta2 <-
  dplyr::filter(df.delta, !is.na(prob) & !is.na(max_pip))


plot_pip_scatter <- function(data,
                             hide.xlab = FALSE,
                             hide.ylab = FALSE,
                             colors = c(
                               "SL" = BuenColors::jdb_palette("brewer_celsius")[8],
                               "NSL" = BuenColors::jdb_palette("brewer_celsius")[2]
                             )) {
  purrr::map(c("SL", "NSL"), function(x) {
    dplyr::mutate(data, color = (prediction == x)) %>%
      dplyr::arrange(color) %>%
      ggplot(aes(max_pip, prob, color = color)) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "grey50"
      ) +
      ggrastr::rasterize(geom_point(), dpi = 300) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25),
        expand = expansion(c(0.05, 0.2))
      ) +
      scale_color_manual(values = c("grey80", colors[x])) +
      locusviz::get_default_theme(hide.xlab = hide.xlab, hide.ylab = (hide.ylab |
        (x == "NSL"))) +
      theme(legend.position = "none") +
      gghighlight::gghighlight(use_direct_label = FALSE) +
      labs(
        x = "Max. PIP (biobanks)",
        y = "PIP (GBMI)",
        title = dplyr::case_when(
          x == "SL" ~ sprintf("%s: Suspicious loci", data$trait[1]),
          x == "NSL" ~ "Non-suspicious loci"
        )
      )
  })
}

plt <-
  dplyr::filter(df.delta2, prob > 0.001 | max_pip > 0.001) %>%
  dplyr::group_split(trait) %>%
  purrr::imap(function(data, i) {
    plts <- plot_pip_scatter(data,
      hide.xlab = ((length(
        unique(df.delta2$trait)
      ) - i) > 1),
      hide.ylab = (i %% 2 == 0)
    )
    plts[[1]] <- plts[[1]] + labs(tag = letters[i])
    return(plts)
  }) %>%
  purrr::flatten() %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = 4)

cowplot::save_plot(
  "figures/SFig8_pip_scatter.pdf",
  plt,
  base_height = 9,
  base_width = 7.2,
  device = cairo_pdf
)