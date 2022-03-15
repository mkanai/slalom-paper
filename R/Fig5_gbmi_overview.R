library(ggplot2)
library(dplyr)
source("~/src/github.com/mkanai/slalom-paper/R/const.R")

df.max_pip_gwas <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/fm_only.max_pip_per_var.b38.tsv.bgz") %>%
  dplyr::select(variant, max_pip_gwas = max_pip)

df.max_pip_eqtl <- rgsutil::read_gsfile("gs://meta-finemapping-simulation/merged.max_pip_per_var.b38.tsv.bgz") %>%
  dplyr::select(variant, max_pip_eqtl = max_pip)


df_dentist <- purrr::map_dfr(gbmi_traits, function(trait) {
  df <- rgsutil::read_gsfile(
    sprintf(
      "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/SLALOM/%s.ABF.snp.bgz",
      trait
    )
  ) %>%
    dplyr::left_join(df.max_pip_gwas) %>%
    dplyr::left_join(df.max_pip_eqtl) %>%
    dplyr::mutate(
      r2 = r**2,
      n_eff_samples = n_samples * (n_cases / n_samples) * (1 - n_cases / n_samples)
    ) %>%
    dplyr::group_by(trait, region) %>%
    dplyr::summarize(
      lead_pip_variant = variant[which.max(prob)],
      n_total = n(),
      n_r2 = sum(r2 > r2_threshold, na.rm = T),
      n_dentist_outlier = sum(
        r2 > r2_threshold &
          nlog10p_dentist > nlog10p_dentist_lower_bound,
        na.rm = T
      ),
      fraction = ifelse(n_r2 == 0, 0, n_dentist_outlier / n_r2),
      n_nonsyn = sum(
        r2 > r2_threshold &
          consequence %in% c("pLoF", "Missense"),
        na.rm = T
      ),
      max_pip = max(prob),
      max_pip_nonsyn = max(prob[r2 > r2_threshold &
        consequence %in% c("pLoF", "Missense")], na.rm = T),
      cs_nonsyn = any(cs[r2 > r2_threshold &
        consequence %in% c("pLoF", "Missense")], na.rm = T),
      cs_99_nonsyn = any(cs_99[r2 > r2_threshold &
        consequence %in% c("pLoF", "Missense")], na.rm = T),
      nonsyn_variants = paste(variant[which(r2 > r2_threshold &
        consequence %in% c("pLoF", "Missense"))], collapse = ","),
      lead_max_pip_gwas = max_pip_gwas[which.max(prob)],
      cs_max_pip_gwas = max(max_pip_gwas[cs], na.rm = T),
      cs_99_max_pip_gwas = max(max_pip_gwas[cs_99], na.rm = T),
      lead_max_pip_eqtl = max_pip_eqtl[which.max(prob)],
      cs_max_pip_eqtl = max(max_pip_eqtl[cs], na.rm = T),
      cs_99_max_pip_eqtl = max(max_pip_eqtl[cs_99], na.rm = T),
      min_neff_r2 = min(n_eff_samples[r2 > r2_threshold], na.rm = T),
      max_neff_r2 = max(n_eff_samples[r2 > r2_threshold], na.rm = T)
    )
  return(df)
})


rgsutil::write_gsfile(
  df_dentist,
  "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi-all-biobank-meta.SLALOM.summary.bgz",
  overwrite = FALSE
)

###################################################################
data <- rgsutil::read_gsfile(
  "gs://meta-finemapping-simulation/gbmi-all-biobank-meta/gbmi-all-biobank-meta.SLALOM.summary.bgz"
) %>%
  dplyr::mutate(
    max_pip_bin = factor(
      cut(max_pip, c(-Inf, 0.1, 0.5, 0.9, 1)),
      levels = c(
        "(-Inf,0.1]",
        "(0.1,0.5]",
        "(0.5,0.9]",
        "(0.9,1]"
      ),
      labels = c("[0,0.1]", "(0.1,0.5]", "(0.5,0.9]", "(0.9,1]"),
    ),
    count_bin = cut(n_dentist_outlier, c(-Inf, 0, 10, 100, Inf)),
    prediction = factor(
      dplyr::case_when(
        max_pip < 0.1 ~ "NA",
        n_dentist_outlier > 0 ~ "SL",
        TRUE ~ "NSL"
      ),
      levels = c("SL", "NSL", "NA")
    ),
    trait = stringr::str_split_fixed(trait, "_", 2)[, 1],
  )

ncol <- 8
max_pip_colors <- c("grey80", BuenColors::jdb_palette("Zissou")[3:5])
names(max_pip_colors) <- c("[0,0.1]", "(0.1,0.5]", "(0.5,0.9]", "(0.9,1]")

traits.order <- dplyr::group_by(data, trait) %>%
  dplyr::summarize(n_loci = length(unique(region))) %>%
  dplyr::arrange(dplyr::desc(n_loci)) %>%
  dplyr::pull(trait)

plt <-
  dplyr::bind_rows(
    dplyr::mutate(data, trait = "All traits"),
    data
  ) %>%
  dplyr::mutate(index = match(trait, c("All traits", traits.order))) %>%
  dplyr::group_by(index) %>%
  dplyr::group_split() %>%
  purrr::map(function(data) {
    i <- data$index[1]
    trait <- data$trait[1]

    df.label <- dplyr::group_by(data, prediction) %>%
      dplyr::summarize(count = n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        frac = count / sum(count),
        label = scales::percent(frac, accuracy = 1)
      )

    ggplot(data, aes(prediction)) +
      geom_bar(aes(fill = max_pip_bin), position = position_stack(reverse = TRUE)) +
      geom_text(
        aes(y = count, label = label),
        data = df.label,
        size = 2,
        vjust = -1
      ) +
      scale_fill_manual(values = max_pip_colors) +
      locusviz::get_default_theme(
        hide.xtitle = TRUE,
        hide.ytitle = (i %% ncol != 1)
      ) +
      locusviz::or_missing(trait != "AAA", theme(legend.position = "none")) +
      labs(title = trait, y = "# loci", fill = "Max PIP bin") +
      scale_x_discrete(drop = FALSE) +
      scale_y_continuous(
        breaks = scales::pretty_breaks(),
        expand = expansion(mult = c(0, 0.2))
      ) +
      theme(plot.title = element_text(
        hjust = ifelse(trait == "Appendectomy", -0.1, 0.1),
        margin = margin(b = -12),
        size = 8
      )) +
      locusviz::or_missing((i %% ncol == 1), theme(plot.tag = element_text(
        face = "bold", margin = margin(r = -6)
      )))
  }) %>%
  append(list(patchwork::guide_area())) %>%
  purrr::reduce(`+`) +
  patchwork::plot_layout(ncol = ncol, guides = "collect") +
  patchwork::plot_annotation(tag_levels = "a")

plt

cowplot::save_plot(
  "figures/Fig5_gbmi_overview.pdf",
  plt,
  base_height = 4,
  base_width = 7.2,
  device = cairo_pdf
)
