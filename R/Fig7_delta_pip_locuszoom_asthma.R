library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

###################################################################
df.max_pip_gwas <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/fm_only.max_pip_per_var.b38.tsv.bgz") %>%
  dplyr::select(variant, max_pip_gwas = max_pip)

df.max_pip_eqtl <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/merged.max_pip_per_var.b38.tsv.bgz") %>%
  dplyr::select(variant, max_pip_eqtl = max_pip)


df.pip <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/bbj_fg_ukbb.shared_trait.pip.tsv.bgz") %>%
  dplyr::mutate(trait = dplyr::case_when(
    trait == "Glaucoma" ~ "POAG",
    trait == "IS" ~ "Stroke",
    TRUE ~ trait
  ))

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

shared_traits <- c(
  "Asthma_Bothsex_inv_var_meta_GBMI_052021",
  "COPD_Bothsex_inv_var_meta_GBMI_052021",
  "Gout_Bothsex_inv_var_meta_GBMI_052021",
  "HF_Bothsex_inv_var_meta_GBMI_052021",
  "IPF_Bothsex_inv_var_meta_GBMI_052021",
  "POAG_Bothsex_inv_var_meta_GBMI_052021",
  "Stroke_Bothsex_inv_var_meta_GBMI_052021",
  "ThC_Bothsex_inv_var_meta_GBMI_052021",
  "VTE_Bothsex_inv_var_meta_GBMI_052021"
)

df <- purrr::map_dfr(shared_traits, function(trait) {
  rgsutil::read_gsfile(
    sprintf(
      "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/SLALOM/%s.ABF.snp.bgz",
      trait
    )
  )
})

df.delta <-
  dplyr::left_join(df, dplyr::select(df.gbmi, trait, region, prediction)) %>%
  dplyr::mutate(trait = stringr::str_split_fixed(trait, "_", 2)[, 1]) %>%
  dplyr::inner_join(df.pip) %>%
  dplyr::left_join(df.max_pip_gwas) %>%
  dplyr::left_join(df.max_pip_eqtl) %>%
  dplyr::mutate(
    ratio = prob / max_pip,
    delta = prob - max_pip,
    prob_bin = cut(prob, pip_bin_breaks),
    max_pip_bin = cut(max_pip, pip_bin_breaks)
  )

rgsutil::write_gsfile(
  df.delta,
  "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi.bbj_fg_ukbb.delta.tsv.bgz",
  overwrite = FALSE
)

##############
df.delta <- rgsutil::read_gsfile(
  "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi.bbj_fg_ukbb.delta.tsv.bgz"
)
df.delta2 <-
  dplyr::filter(df.delta, !is.na(prob) & !is.na(max_pip))

count_data <- function(data, threshold = 0.99) {
  dplyr::bind_rows(
    dplyr::arrange(data, prob) %>%
      dplyr::mutate(
        # q = prob,
        q = dplyr::row_number() / n(),
        q_bin = threshold,
        cohort = "GBMI"
      ) %>%
      dplyr::filter(q > threshold),
    dplyr::arrange(data, max_pip) %>%
      dplyr::mutate(
        # q = max_pip,
        q = dplyr::row_number() / n(),
        q_bin = threshold,
        cohort = "biobanks"
      ) %>%
      dplyr::filter(q > threshold),
  ) %>%
    dplyr::mutate(
      consequence = dplyr::case_when(
        consequence %in% c("pLoF", "Missense", "Synonymous") ~ "Coding",
        consequence %in% c("UTR5", "UTR3") ~ "UTR",
        is.na(consequence) ~ NA_character_,
        TRUE ~ consequence
      ),
      consequence = factor(
        consequence,
        levels = c(
          "Coding",
          "UTR",
          "Promoter",
          "CRE",
          "Non-genic"
        )
      )
    ) %>%
    dplyr::filter(!is.na(consequence) & !is.na(q_bin)) %>%
    dplyr::group_by(cohort, q_bin, consequence) %>%
    dplyr::count() %>%
    dplyr::group_by(cohort, q_bin) %>%
    dplyr::mutate(
      total = sum(n),
      n2 = total - n
    ) %>%
    dplyr::ungroup()
}

plot_enrichment_ranking <- function(data,
                                    title,
                                    hide.ylab = FALSE,
                                    ylim = c(0.8, 5)) {
  pd <- position_dodge(width = 0.75)
  purrr::map_dfr(
    c(0.995, 0.999, 0.9995),
    function(x) {
      count_data(data, threshold = x)
    }
  ) %>%
    dplyr::group_split(q_bin, consequence) %>%
    purrr::map_dfr(function(data) {
      print(data)
      if (nrow(data) == 1) {
        return(NULL)
      }

      m <-
        magrittr::extract(data, match(c("biobanks", "GBMI"), data$cohort), c("n2", "n")) %>%
        as.matrix() %>%
        epitools::riskratio(method = "boot")
      print(m)

      tibble::tibble(
        q_bin = data$q_bin[1],
        consequence = data$consequence[1],
        enrichment = m$measure[2, "estimate"],
        lower = m$measure[2, "lower"],
        upper = m$measure[2, "upper"]
      )
    }) %>%
    dplyr::mutate(
      q_bin = factor(q_bin),
      q_bin = forcats::fct_recode(
        q_bin,
        "Top 5%" = "0.95",
        "Top 1%" = "0.99",
        "Top 0.5%" = "0.995",
        "Top 0.1%" = "0.999",
        "Top 0.05%" = "0.9995",
        "Top 0.01%" = "0.9999",
        "Top 0.005%" = "0.99995",
      )
    ) %>%
    ggplot(aes(consequence,
      enrichment,
      color = consequence,
      # color = q_bin
      shape = q_bin
    )) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_point(position = pd) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
      width = 0,
      position = pd
    ) +
    labs(color = "PIP ranking threshold") +
    locusviz::get_default_theme(
      hide.xtitle = TRUE,
      hide.ylab = hide.ylab
    ) +
    scale_color_manual(values = unname(annot_colors[c("pLoF", "UTR5", "Promoter", "CRE", "Non-genic")]), guide = FALSE) +
    coord_cartesian(ylim = ylim) +
    labs(y = "Enrichment", title = title, shape = "PIP ranking")
}

plts <-
  list(
    plot_enrichment_ranking(df.delta2, title = "All loci") +
      labs(tag = "a"),
    dplyr::filter(df.delta2, prediction == "SL") %>%
      plot_enrichment_ranking(
        hide.ylab = TRUE,
        title = "Suspicious loci"
      ) +
      theme(legend.position = "none") +
      labs(tag = "b"),
    dplyr::filter(df.delta2, prediction == "NSL") %>%
      plot_enrichment_ranking(
        hide.ylab = TRUE,
        title = "Non-suspicious loci"
      ) +
      theme(legend.position = "none") +
      labs(tag = "c")
  )

###################################################################
dplyr::summarize(
  df.delta2,
  n_cs_nonsyn = sum(cs &
    (consequence %in% c("pLoF", "Missense"))),
  n_cs = sum(cs),
  n_cs_any_nonsyn = sum(any_cs &
    (consequence %in% c("pLoF", "Missense"))),
  n_cs_any = sum(any_cs)
)

###################################################################

fm_insights <-
  reticulate::import_from_path("fm_insights", path = "~/src/github.com/mkanai/finemapping-insights")

df.csm <-
  rgsutil::read_gsfile(fm_insights$get_merged_results_path("fm_only.csm_id", "tsv.bgz")) %>%
  dplyr::select(cohort, trait, region, variant, susie.cs_id, csm_id)

df.position_b38 <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/fm_only.max_pip_per_var.b38.tsv.bgz") %>%
  dplyr::mutate(
    position = locusviz::parse_variant(variant)$position,
    variant = variant_b37
  ) %>%
  dplyr::select(variant, position)


plot_multi_locuszoom_with_gbmi <- function(df.gbmi,
                                           trait,
                                           lead_variant,
                                           lead_variant_b37,
                                           rsid,
                                           window,
                                           trait_label = NULL,
                                           cohorts = c("GBMI", "BBJ", "FG", "UKBB"),
                                           df.csm = NULL,
                                           additional_highlight_pos = NULL,
                                           additional_highlight_labels = NULL,
                                           hide.trait = FALSE,
                                           fm.height = 0.3,
                                           fontsize = 8) {
  lead_pos <- locusviz::parse_variant(lead_variant)$position
  chromosome <- locusviz::parse_variant(lead_variant)$chromosome
  lead_pos_b37 <- locusviz::parse_variant(lead_variant_b37)$position
  chromosome_b38 <-
    locusviz::parse_variant(lead_variant_b37)$chromosome

  if (length(window) == 1) {
    window <- rep(window, 2)
  }
  start <- lead_pos - window[1]
  end <- lead_pos + window[2]
  start_b37 <- lead_pos_b37 - window[1]
  end_b37 <- lead_pos_b37 + window[2]

  highlight_pos <- lead_pos
  if (!is.null(additional_highlight_pos)) {
    highlight_pos <- c(highlight_pos, additional_highlight_pos)
  }

  if (is.null(trait_label)) {
    trait_label <- trait
  }

  lead_csm_id <-
    dplyr::filter(df.csm, trait == .env$trait &
      variant == lead_variant_b37) %>%
    dplyr::pull(csm_id) %>%
    head(1)

  # dirty hack for now
  if (trait == "Asthma" & lead_variant_b37 == "1:152285861:G:A") {
    lead_csm_id <- 4
  }

  if (length(lead_csm_id) == 1) {
    df.gbmi <- dplyr::mutate(df.gbmi, cs = ifelse(cs == 1, lead_csm_id, cs))
  } else {
    lead_csm_id <- 1
  }

  cs.colors <-
    dplyr::filter(df.csm, trait == .env$trait &
      cohort %in% cohorts) %>%
    dplyr::bind_rows(
      tibble::tibble(
        cohort = "GBMI",
        trait = .env$trait,
        variant = lead_variant_b37,
        csm_id = lead_csm_id
      )
    ) %>%
    dplyr::mutate(locusviz::parse_variant(variant)) %>%
    dplyr::filter(
      chromosome == .env$chromosome_b38 &
        .env$start_b37 * 0.99 <= position &
        position <= .env$end_b37 * 1.01
    ) %>%
    dplyr::arrange(abs(position - position[which(.$variant == lead_variant_b37)[1]])) %$%
    locusviz::get_cs_color_mapping(.$csm_id, highlight_cs_ids = .$csm_id[which(.$variant == lead_variant_b37)])

  print(cs.colors)
  print(table(df.gbmi$cs))


  panels <-
    purrr::flatten(purrr::map(cohorts, function(cohort) {
      if (cohort == "GBMI") {
        data <-
          locusviz::preprocess(df.gbmi, pip_col = "prob", cs_id_col = "cs")
      } else {
        data <-
          rgsutil::read_gsfile(fm_insights$get_locus_txt_path(cohort, trait, lead_variant_b37)) %>%
          dplyr::left_join(df.csm) %>%
          dplyr::mutate(locusviz::parse_variant(variant)) %>%
          dplyr::select(-position) %>%
          dplyr::left_join(df.position_b38) %>%
          locusviz::preprocess(
            lead_variant = lead_variant_b37,
            beta_col = "beta_marginal",
            se_col = "se_marginal",
            cs_id_col = "csm_id",
            r2_col = "lead_r2"
          )
      }

      p_manhattan <- locusviz::plot_manhattan_panel(
        data,
        highlight_pos = highlight_pos,
        xlim = c(start, end),
        plot.loglog_p = FALSE,
        nlog10p_threshold = 0,
        point.size = 0.5,
        point.size2 = 3,
        line.size = 0.2,
        title = ifelse(hide.trait, cohort, paste(cohort, trait_label, sep = ": ")),
        ggtheme = locusviz::get_default_theme(
          fontsize = fontsize,
          hide.xtext = TRUE,
          hide.xtitle = TRUE,
          legend.position = if (cohort == cohorts[1]) {
            c(1, 1.1)
          } else {
            "none"
          }
        ) + theme(
          plot.title = element_text(
            hjust = 0.01,
            margin = margin(b = -12),
            size = fontsize
          ),
          legend.margin = margin(0, 0, 0, 0),
          # legend.title = element_text(size = fontsize),
          # legend.key.size = unit(0.15, units = "cm")
        ),
        rasterize = TRUE
      )

      if (cohort == cohorts[1]) {
        label_pos <- data$nlog10p[which(data$position == lead_pos)[1]]
        if (label_pos < 10) {
          label_pos <- label_pos + max(data$nlog10p, na.rm = T) / 10
        }
        p_manhattan <- p_manhattan +
          geom_text(
            aes(x, y, label = label),
            data = tibble::tibble(
              x = lead_pos,
              y = label_pos,
              label = rsid
            ),
            size = 2,
            vjust = -1
          ) + scale_y_continuous(expand = expansion(c(0, 0.2)))
        if (!is.null(additional_highlight_labels)) {
          p_manhattan <- p_manhattan +
            purrr::pmap(tibble::tibble(position = additional_highlight_pos, label = additional_highlight_labels), function(position, label) {
              label_pos <- data$nlog10p[which(data$position == position)[1]]
              if (label_pos < 10) {
                label_pos <- label_pos + max(data$nlog10p, na.rm = T) / 10
              }
              geom_text(
                aes(x, y, label = label),
                data = tibble::tibble(
                  x = .env$position,
                  y = label_pos,
                  label = .env$label
                ),
                size = 2,
                vjust = -1
              )
            })
        }
      }

      p_fm <- locusviz::plot_fm_panel(
        data,
        highlight_pos = highlight_pos,
        xlim = c(start, end),
        ylim = c(0, 1.1),
        ybreaks = seq(0, 1, length = 2),
        point.size = 0.5,
        point.size2 = 3,
        ggtheme = locusviz::get_default_theme(
          fontsize = fontsize,
          hide.xtext = cohort != cohorts[length(cohorts)],
          hide.xtitle = TRUE,
          legend.position = "none"
        ) + theme(plot.margin = margin(0, 0.1, 0.2, 0.1, unit = "cm")),
        rasterize = TRUE,
        cs.colors = cs.colors
      )

      return(list(p_manhattan, p_fm))
    }))

  panels <- c(panels, list(
    locusviz::plot_gene_panel(
      chromosome,
      start,
      end,
      genome_build = "hg38",
      highlight_pos = highlight_pos,
      fontsize = fontsize,
      point.size = 1.5,
      label.size = 1.5,
      length = unit(0.05, "cm")
    )
  ))

  plt <-
    purrr::reduce(c(panels), `+`) + patchwork::plot_layout(ncol = 1, heights = c(rep(c(1, fm.height), length(cohorts)), 0.1))
  return(plt)
}

# IL33 improvement
p1 <- plot_multi_locuszoom_with_gbmi(
  df.gbmi = munge_gbmi_sumstats(gbmi_traits["Asthma"], "chr9:3629917-7270078"),
  trait = "Asthma",
  lead_variant = "chr9:6197392:T:C",
  lead_variant_b37 = "9:6197392:T:C",
  rsid = "rs1888909",
  window = 60000,
  cohort = c("GBMI", "FG"),
  hide.trait = TRUE,
  df.csm = df.csm
)

# OTULINL missense
p2 <- plot_multi_locuszoom_with_gbmi(
  df.gbmi = munge_gbmi_sumstats(gbmi_traits["Asthma"], "chr5:14110200-15110200"),
  trait = "Asthma",
  lead_variant = "chr5:14610200:C:G",
  lead_variant_b37 = "5:14610309:C:G",
  rsid = "rs16903574",
  window = 50000,
  cohort = c("GBMI", "UKBB"),
  additional_highlight_pos = c(14633803),
  additional_highlight_labels = c("rs528167451"),
  hide.trait = TRUE,
  df.csm = df.csm
)

# IL13 improvement
p3 <- plot_multi_locuszoom_with_gbmi(
  df.gbmi = munge_gbmi_sumstats(gbmi_traits["Asthma"], "chr5:130224057-133160151"),
  trait = "Asthma",
  lead_variant = "chr5:132660151:T:C",
  lead_variant_b37 = "5:131995843:T:C",
  rsid = "rs1295686",
  window = c(2000, 1200),
  additional_highlight_pos = c(132660272),
  additional_highlight_labels = c(""),
  cohort = c("GBMI", "UKBB"),
  hide.trait = TRUE,
  df.csm = df.csm
)
p3[[1]] <- p3[[1]] + geom_text(
  aes(x, y, label = label),
  data = tibble::tibble(
    x = 132660272,
    y = -log10(2.442e-47),
    label = "rs20541"
  ),
  size = 2,
  hjust = -0.2
)

# FLG, multi causal
# TODO: export rs12123821, not current FLG stop-gained
p4 <- plot_multi_locuszoom_with_gbmi(
  df.gbmi = munge_gbmi_sumstats(gbmi_traits["Asthma"], "chr1:151706676-153299487"),
  trait = "Asthma",
  lead_variant = "chr1:152313385:G:A",
  lead_variant_b37 = "1:152285861:G:A",
  rsid = "rs61816761",
  window = c(175000, 50000),
  cohort = c("GBMI", "UKBB"),
  additional_highlight_pos = c(152206676),
  additional_highlight_labels = c("rs12123821"),
  hide.trait = TRUE,
  df.csm = df.csm
)



p1[[1]] <- p1[[1]] + labs(tag = "d")
p2[[1]] <- p2[[1]] + labs(tag = "e")
p3[[1]] <- p3[[1]] + labs(tag = "f")
p4[[1]] <- p4[[1]] + labs(tag = "g")

plt <- purrr::reduce(plts, `|`) / (p1 |
  p2) / (p3 |
  p4) + patchwork::plot_layout(heights = c(2, 3.0, 3.0))

cowplot::save_plot(
  "figures/Fig7_delta_pip_locuszoom_asthma.pdf",
  plt,
  base_height = 8,
  base_width = 7.2,
  device = cairo_pdf
)
