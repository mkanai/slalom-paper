nlog10p_dentist_lower_bound <- 4
r2_threshold <- 0.6
count_threshold <- 1
pip_bin_breaks <- c(0, 0.01, 0.1, 0.5, 0.9, 1.0)

gnomad_pop_colors <- c(
  AFR = "#941494",
  AMR = "#ED1E24",
  EAS = "#108C44",
  FIN = "#002F6C",
  NFE = "#6AA5CD",
  SAS = "#FF9912",
  Meta = "black"
)

gbmi_traits <- c(
  "AAA" = "AAA_Bothsex_inv_var_meta_GBMI_052021",
  "Asthma" = "Asthma_Bothsex_inv_var_meta_GBMI_052021",
  "AcApp" = "AcApp_Bothsex_inv_var_meta_GBMI_052021",
  "Appendectomy" = "Appendectomy_Bothsex_inv_var_meta_GBMI_052021",
  "COPD" = "COPD_Bothsex_inv_var_meta_GBMI_052021",
  "Gout" = "Gout_Bothsex_inv_var_meta_GBMI_052021",
  "HCM" = "HCM_Bothsex_inv_var_meta_GBMI_052021",
  "HF" = "HF_Bothsex_inv_var_meta_GBMI_052021",
  "IPF" =  "IPF_Bothsex_inv_var_meta_GBMI_052021",
  "POAG" = "POAG_Bothsex_inv_var_meta_GBMI_052021",
  "Stroke" =  "Stroke_Bothsex_inv_var_meta_GBMI_052021",
  "ThC" = "ThC_Bothsex_inv_var_meta_GBMI_052021",
  "UtC" = "UtC_Female_inv_var_meta_GBMI_052021",
  "VTE" = "VTE_Bothsex_inv_var_meta_GBMI_052021"
)

# functional annotations
annot_levels <- c(
  "pLoF",
  "Missense",
  "Synonymous",
  "UTR5",
  "UTR3",
  "Promoter",
  "CRE",
  "Non-genic"
)
annot_colors <- c(
  BuenColors::jdb_palette("brewer_celsius")[c(8, 6)],
  "#AAAAAA",
  BuenColors::jdb_palette("brewer_purple")[c(6, 4)],
  BuenColors::jdb_palette("brewer_blue")[c(8, 6, 3)]
)
names(annot_colors) <- annot_levels

na_binconf <- function(x, n, lower_upper, ...) {
  ret <- rep(NA, length(x))
  idx <- !is.na(x) & !is.na(n)
  if (sum(idx) > 0) {
    ret[idx] <- Hmisc::binconf(x[idx], n[idx], return.df = TRUE, ...)[[lower_upper]]
  }
  return(ret)
}

munge_gbmi_sumstats <- function(trait, region) {
  data <- rgsutil::read_gsfile(
    sprintf(
      "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/SLALOM/%s.ABF.snp.bgz",
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

plot_locuszoom_slalom <- function(data,
                                  trait,
                                  region,
                                  title = NULL,
                                  additional_lead_pos = NULL,
                                  window_size = 50000,
                                  fontsize = 6) {
  data <- locusviz::preprocess(data, pip_col = "prob", cs_id_col = "cs")
  region <- data$region[1]
  lead_idx <- which.max(data$pip)
  lead_pos <- data$position[lead_idx]
  lead_variant <- data$variant[lead_idx]
  lead_rsid <- data$rsid[lead_idx]
  lead_z <- data$beta[lead_idx] / data$se[lead_idx]
  lead_nlog10p <- data$nlog10p[lead_idx]
  lead_n_samples <- data$n_samples[lead_idx]
  chromosome <- data$chromosome[1]

  ld_range <- dplyr::filter(data, r2 > 0.6) %>%
    dplyr::pull(position) %>%
    range() %>%
    diff()

  window <- min(max(ld_range + 25000, 50000), window_size)
  start <- lead_pos - window
  end <- lead_pos + window
  p_test_threshold <- 10**-nlog10p_dentist_lower_bound

  if (min(data$p) > 5e-8 |
    length(lead_idx) == 0) {
    print(region)
    return(NULL)
  }

  if (is.null(title)) {
    title <- sprintf(
      "%s",
      stringr::str_split_fixed(trait, "_", 2)[, 1]
    )
  }

  n_outlier_variants <- with(data, sum(outlier,
    na.rm = T
  ))

  p1 <-
    locusviz::plot_locuszoom(
      data,
      highlight_pos = c(lead_pos, additional_lead_pos),
      xlim = c(start, end),
      gene.args = list(genome_build = "hg38"),
      manhattan.title = title,
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
      fontsize = fontsize,
      rasterize = TRUE,
      patchwork = FALSE
    )
  p1[[1]] <- p1[[1]] + geom_text(
    aes(
      x = lead_pos,
      y = lead_nlog10p,
      label = lead_rsid
    ),
    data = tibble::tibble(),
    size = 2,
    vjust = -1
  )
  p1[[3]] <- p1[[3]] + theme(legend.margin = margin())

  p2 <-
    ggplot() +
    geom_hline(
      yintercept = qchisq(5e-8, 1, lower.tail = F),
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_abline(
      intercept = 0,
      slope = data$chisq[lead_idx],
      linetype = "dashed",
      color = "grey50"
    ) +
    ggrastr::rasterize(geom_point(
      aes(
        r2,
        chisq,
        color = nlog10p_dentist,
        shape = outlier,
        size = outlier
      ),
      data = dplyr::arrange(data, nlog10p_dentist) %>% dplyr::filter(!outlier)
    ),
    dpi = 300
    ) +
    geom_point(
      aes(
        r2,
        chisq,
        color = nlog10p_dentist,
        shape = outlier,
        size = outlier
      ),
      data = dplyr::filter(data, outlier)
    ) +
    geom_point(
      aes(x, y),
      color = "purple3",
      shape = 18,
      size = 3,
      data = tibble::tibble(x = 1, y = lead_z**2)
    ) +
    scale_color_gradientn(
      colors = BuenColors::jdb_palette("solar_extra"),
      limits = c(0, nlog10p_dentist_lower_bound),
      oob = scales::squish
    ) +
    scale_shape_manual(
      values = c(`TRUE` = 18, `FALSE` = 20),
      guide = FALSE
    ) +
    scale_size_manual(
      values = c(`TRUE` = 3, `FALSE` = 1.5),
      guide = FALSE
    ) +
    locusviz::get_default_theme(
      fontsize = fontsize,
      legend.position = c(1, 0),
      legend.justification = c(1, 0)
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.02,
        margin = margin(b = -12)
      ),
      legend.direction = "horizontal",
      legend.key.width = unit(0.2, "cm")
    ) +
    labs(
      x = expression(paste(italic(r)^2, "to the lead variant")),
      y = expression(italic(chi)^2),
      color = expression(paste(-log[10], " ", italic(P)[`DENTIST-S`])),
      title = sprintf(
        "%s DENTIST-S outlier variant%s",
        ifelse(n_outlier_variants > 0, n_outlier_variants, "No"),
        ifelse(n_outlier_variants > 1, "s", "")
      )
    ) +
    scale_x_continuous(limits = c(0, 1))

  return(list(p1, p2))
}

##############################################################################
# Fig. 4, 6, SFig. 6
##############################################################################
plot_validation_fraction <- function(data,
                                     title,
                                     ylab = NULL,
                                     hide.ylab = FALSE,
                                     hide.legend = FALSE) {
  pd <- position_dodge(width = 0.75)

  data2 <-
    dplyr::filter(data, prediction != "NA") %>%
    dplyr::mutate(prediction = factor(
      prediction,
      levels = c("SL", "NSL"),
      labels = c("Suspicious loci", "Non-suspicious loci")
    )) %>%
    dplyr::group_by(prediction) %>%
    dplyr::summarize(
      count = colSums(cbind(count_lead, count_cs, count_cs_99), na.rm = T),
      count2 = colSums(cbind(!count_lead, !count_cs, !count_cs_99), na.rm = T),
      type = factor(
        c("Lead PIP", "in 95% CS", "in 99% CS"),
        levels = c("Lead PIP", "in 95% CS", "in 99% CS")
      ),
      frac = count / (count + count2),
      frac.lower = na_binconf(count, rep(n(), 3), "Lower"),
      frac.upper = na_binconf(count, rep(n(), 3), "Upper")
    ) %>%
    dplyr::ungroup()

  g_bracket <-
    dplyr::group_split(data2, type) %>%
    purrr::map_dfr(function(x) {
      m <- dplyr::select(x, count2, count) %>%
        as.matrix() %>%
        magrittr::extract(1:2, 1:2) %>%
        epitools::riskratio(method = "boot")
      print(m)
      pvalue <- m$p.value[2, 2]
      rr <- m$measure[2, "estimate"]
      tip.length <- 0.05
      scale <- (max(data2$frac.upper) - min(data2$frac.lower))
      return(
        tibble::tibble(
          type = x$type[1],
          pvalue = pvalue,
          rr = rr,
          y.position = max(x$frac.upper) + tip.length,
          tip.length1 = (y.position - x$frac.upper[1] - 0.04) / scale,
          tip.length2 = (y.position - x$frac.upper[2] - 0.04) / scale
        )
      )
    }) %>%
    dplyr::filter(pvalue < 0.05) %>%
    purrr::pmap(function(type,
                         pvalue,
                         rr,
                         y.position,
                         tip.length1,
                         tip.length2) {
      index <- as.numeric(type) - 1
      list(
        ggpubr::geom_bracket(
          xmin = index + 0.81,
          xmax = index + 1.19,
          y.position = y.position,
          # tip.length = c(tip.length1, tip.length2),
          color = "black",
          label = symnum(
            pvalue,
            cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
            symbols = c("****", "***", "**", "*", "ns")
          )
        ),
        geom_text(
          aes(
            x = index + 1,
            y = y.position,
            label = sprintf(
              "%sx",
              scales::number(rr, accuracy = 0.1)
            )
          ),
          color = "black",
          data = tibble::tibble(),
          vjust = -3,
          size = 2
        )
      )
    })

  ggplot(data2, aes(type, frac, color = prediction)) +
    geom_point(position = pd) +
    geom_errorbar(aes(ymin = frac.lower, ymax = frac.upper),
      width = 0,
      position = pd
    ) +
    g_bracket +
    locusviz::get_default_theme(hide.xtitle = TRUE, hide.ylab = hide.ylab) +
    theme(
      legend.position = c(1, 0.1),
      legend.justification = c(1, 0),
      plot.title = element_text(
        hjust = 0.05,
        margin = margin(b = -12 * (
          stringr::str_count(title, "\n") + 1
        ))
      )
    ) +
    locusviz::or_missing(hide.legend, theme(legend.position = "none")) +
    scale_color_manual(values = BuenColors::jdb_palette("brewer_celsius")[c(8, 2)]) +
    labs(y = ylab, color = "Prediction", title = title) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      expand = expansion(c(0, 0.2), 0),
      labels = scales::percent_format()
    ) +
    coord_cartesian(ylim = c(0, 1))
}

plot_validataion_wrap <- function(data, cohort) {
  plts <- list(
    dplyr::mutate(
      data,
      count_lead = max_pip == max_pip_nonsyn,
      # count_lead = max_pip_nonsyn > 0.01 & max_pip != max_pip_nonsyn,
      count_cs = cs_nonsyn,
      count_cs_99 = cs_99_nonsyn
    ) %>%
      dplyr::filter(n_nonsyn > 0) %>%
      plot_validation_fraction(
        title = "Nonsynonymous coding variants",
        ylab = sprintf("%% %s locus\n(with tagging nonsyn. & max PIP > 0.1)", cohort)
      ),
    dplyr::mutate(
      data,
      count_lead = lead_max_pip_gwas > 0.9,
      count_cs = cs_max_pip_gwas > 0.9,
      count_cs_99 = cs_99_max_pip_gwas > 0.9
    ) %>%
      plot_validation_fraction(
        title = "High-PIP (> 0.9) complex trait variants\n                       in biobank fine-mapping",
        ylab = sprintf("%% %s locus (max PIP > 0.1)", cohort),
        hide.legend = TRUE
      ),
    dplyr::mutate(
      data,
      count_lead = lead_max_pip_eqtl > 0.9,
      count_cs = cs_max_pip_eqtl > 0.9,
      count_cs_99 = cs_99_max_pip_eqtl > 0.9
    ) %>%
      plot_validation_fraction(
        title = "High-PIP (> 0.9) cis-eQTL variants\n            in GTEx and eQTL catalogue",
        hide.ylab = TRUE,
        hide.legend = TRUE
      )
  )
  return(plts)
}

##############################################################################
# SFig. 5,7
##############################################################################
plot_sample_size_ratio_distribution <- function(data,
                                                fill,
                                                title = NULL,
                                                label.y = 0.8,
                                                hjust = -0.1,
                                                xlim = NULL,
                                                log10 = FALSE) {
  m <- median(data$ratio, na.rm = T)
  if (log10) {
    scale_x <- scale_x_log10(limits = xlim, expand = expansion(0, 0))
  } else {
    scale_x <- scale_x_continuous(limits = xlim, expand = expansion(0, 0))
  }
  ggplot(data, aes(ratio)) +
    ggridges::geom_density_line(fill = fill) +
    geom_vline(
      xintercept = m,
      color = "grey50",
      linetype = "dashed"
    ) +
    geom_text(
      aes(
        x = m,
        y = label.y,
        label = sprintf("median = %s", scales::number(m, accuracy = 0.1))
      ),
      size = 2,
      hjust = hjust,
      data = tibble::tibble()
    ) +
    scale_x +
    scale_y_continuous(expand = expansion(c(0, 0.1), 0)) +
    locusviz::get_default_theme() +
    labs(x = "Sample size ratio", y = "Density", title = title)
}