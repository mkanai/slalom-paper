library(dplyr)
library(ggplot2)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

munge <- function(trait, region) {
  data <- rgsutil::read_gsfile(
    sprintf(
      "gs://meta-finemapping-simulation/covid19_hgi/results/%s.ABF.snp.bgz",
      trait
    )
  ) %>%
    dplyr::filter(region == .env$region) %>%
    dplyr::mutate(
      chisq = (beta / se)**2,
      r2 = r**2,
      cs = ifelse(cs, 1, -1),
      outlier = r2 > r2_threshold &
        nlog10p_dentist > nlog10p_dentist_lower_bound
    )
}

# release 5
trait.r5 <- "COVID19_HGI_B2_ALL_20210107"
lead_variant.r5 <- "chr19:10317045:T:A"
lead_pos.r5 <- 10317045
region.r5 <- "chr19:9817045-10817045"


# rs34536443
lead_variant2 <- "chr19:10352442:G:C"
lead_pos2 <- 10352442

chromosome <- "chr19"
start <- lead_pos2 - 100000
end <- lead_pos2 + 50000

data.r5 <- munge(trait.r5, region.r5) %>%
  locusviz::preprocess(pip_col = "prob", cs_id_col = "cs")

lead_idx <- which(data.r5$variant %in% c(lead_variant.r5, lead_variant2))

p1 <-
  locusviz::plot_locuszoom(
    data.r5,
    highlight_pos = c(lead_pos.r5, lead_pos2),
    xlim = c(start, end),
    gene.args = list(genome_build = "hg38"),
    manhattan.title = "COVID-19 hospitalization (release 5)",
    manhattan.loglog_p = FALSE,
    nlog10p_threshold = 0,
    manhattan.args = list(
      point.size = 0.5,
      point.size2 = 3,
      line.size = 0.2
    ),
    fm.args = list(
      point.size = 0.5,
      point.size2 = 3
    ),
    r2.args = list(
      ybreaks = seq(0, 1, by = 0.5),
      point.size = 0.5,
      point.size2 = 3
    ),
    fm.ylim = c(0, 1),
    fm.breaks = seq(0, 1, by = 0.5),
    plot.fm = TRUE,
    plot.r2 = TRUE,
    plot.gene_score = FALSE,
    plot.gene = TRUE,
    fontsize = 8,
    rasterize = TRUE,
    patchwork = FALSE
  )

p1[[1]] <- p1[[1]] +
  geom_text(
    aes(
      x = data.r5$position[lead_idx],
      y = -log10(data.r5$p[lead_idx]),
      label = c("rs74956615", data.r5$rsid[lead_idx[2]])
    ),
    data = tibble::tibble(),
    size = 2,
    vjust = -1
  )

plt <- purrr::reduce(p1, `+`) + patchwork::plot_layout(
  heights = c(1, 0.3, 0.3, 0.05),
  ncol = 1
)
plt

cowplot::save_plot(
  "figures/SFig1_COVID19_HGI_TYK2_locuszoom.pdf",
  plt,
  base_height = 4,
  base_width = 7.2,
  device = cairo_pdf
)
