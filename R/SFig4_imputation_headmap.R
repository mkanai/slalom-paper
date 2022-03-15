library(ggplot2)
library(dplyr)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

df.shared.grid <- purrr::pmap_dfr(
  tidyr::crossing(
    maf_str = c("0001", "0005", "001", "005"),
    rsq_str = c("02", "04", "06", "08")
  ),
  function(maf_str, rsq_str) {
    df.count <- rgsutil::read_gsfile(
      sprintf(
        "gs://meta-finemapping-simulation/imputed_info/upset.gnomad_min_maf%s.rsq%s.tsv.bgz",
        maf_str,
        rsq_str
      )
    )

    dplyr::summarize(
      df.count,
      maf_str = maf_str,
      rsq_str = rsq_str,
      numer = n[which(collection == "abcdefghijklmnopqrstu")],
      denom = sum(n),
      frac = numer / denom
    )
  }
) %>%
  dplyr::mutate(
    maf = factor(stringr::str_replace(maf_str, "^0", "0.")),
    rsq = factor(stringr::str_replace(rsq_str, "^0", "0."))
  )

plt <-
  list(
    ggplot(df.shared.grid, aes(maf, rsq)) +
      geom_tile(aes(fill = frac)) +
      geom_text(aes(label = scales::percent(frac, accuracy = 0.1), color = rsq), size = 2) +
      scale_fill_viridis_c(label = scales::percent_format(accuracy = 1)) +
      scale_color_manual(values = c(rep("black", 3), "white"), guide = FALSE) +
      labs(
        x = "gnomAD MAF threshold",
        y = "Rsq threshold",
        fill = "% shared",
        title = "Fraction of shared variants in every combination"
      ) +
      locusviz::get_default_theme(legend.position = "bottom") +
      scale_x_discrete(expand = expansion()) +
      scale_y_discrete(expand = expansion(c(0, 0.25))) +
      theme(
        legend.key.width = unit(1.2, "cm"),
        plot.title = element_text(hjust = 0.02)
      ),
    ggplot(df.shared.grid, aes(maf, rsq)) +
      geom_tile(aes(fill = numer)) +
      geom_text(aes(label = scales::label_number_si()(numer), color = maf), size = 2) +
      scale_fill_viridis_c(option = "magma") +
      scale_color_manual(values = c(rep("black", 3), "white"), guide = FALSE) +
      labs(
        x = "gnomAD MAF threshold",
        y = "Rsq threshold",
        fill = "# shared",
        title = "No. shared variants in every combination"
      ) +
      locusviz::get_default_theme(legend.position = "bottom", hide.ylab = TRUE) +
      scale_x_discrete(expand = expansion()) +
      scale_y_discrete(expand = expansion(c(0, 0.25))) +
      theme(
        legend.key.width = unit(1.2, "cm"),
        plot.title = element_text(hjust = 0.02)
      )
  ) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(tag_levels = "a")

plt

cowplot::save_plot(
  "figures/SFig4_shared_heatmap.pdf",
  plt,
  base_height = 4.2,
  base_width = 7.2,
  device = cairo_pdf
)
