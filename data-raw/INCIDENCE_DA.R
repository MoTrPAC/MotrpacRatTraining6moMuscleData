library(MotrpacRatTraining6moMuscleData)
library(dplyr)
library(data.table)


# Nested list of DA results
DA_list <- list(
  "TRNSCRPT_GN" = TRNSCRPT_GN_DA,
  "PROT_GN" = PROT_GN_DA,
  "PHOSPHO_GN" = PHOSPHO_GN_NORM_DA,
  "ACETYL_GN" = ACETYL_GN_NORM_DA,
  "REDOX_GN" = REDOX_GN_NORM_DA,
  "METAB_GN" = METAB_GN_DA,
  "TRNSCRPT_VL" = TRNSCRPT_VL_DA,
  "METAB_VL" = METAB_VL_DA
) %>%
  purrr::transpose()

INCIDENCE_DA <- lapply(DA_list, function(group) {
  res <- lapply(group, function(ome) {
    out <- ome %>%
      filter(adj.P.Val < 0.05) %>%
      select(
        contrast, featureName,
        any_of(c(
          "gene_symbol" = "gene_symbol",
          "gene_symbol" = "external_gene_name"
        ))
      )

    if (!"gene_symbol" %in% colnames(out)) {
      out$gene_symbol <- NA
    }

    return(out)
  }) %>%
    bind_rows(.id = "ome") %>%
    mutate(
      tissue = sub(".*_", "", ome),
      tissue = factor(tissue, levels = c("GN", "VL")),
      ome = sub("_.*", "", ome),
      # Ordered according to central dogma
      ome = factor(ome, levels = c(
        "TRNSCRPT", "PROT", "PHOSPHO",
        "ACETYL", "REDOX", "METAB"
      )),
      is_signif = 1L
    ) %>%
    tidyr::pivot_wider(
      id_cols = any_of(c("tissue", "ome", "featureName", "gene_symbol")),
      names_from = contrast,
      values_from = is_signif,
      values_fill = 0L
    )

  if (sum(grepl(" - ", colnames(res))) > 2L) {
    res <- res %>%
      # Keys for searching specific intersections quickly
      mutate(
        key_female = apply(
          across(.cols = starts_with("F_")), 1,
          \(.i) paste(c("F", .i), collapse = "")
        ),
        key_male = apply(
          across(.cols = starts_with("M_")), 1,
          \(.i) paste(c("M", .i), collapse = "")
        )
      )
  } else {
    res <- res %>%
      mutate(key = apply(
        across(.cols = starts_with("M_")), 1,
        \(.i) paste(c("X", .i), collapse = "")
      ))
  }

  res <- res %>%
    mutate(across(
      .cols = contains("key"),
      .fns = ~ factor(.x, levels = sort(unique(.x)))
    )) %>%
    # Move contrast columns to end
    relocate(contains("-"), .after = everything())

  # Convert to data.table to reorder rows if some columns are not present
  setDT(res)

  col_order <- c(
    "tissue" = 1, "ome" = 1, "key_female" = 1, "key_male" = 1,
    "key" = 1, "featureName" = 1, "gene_symbol" = 1
  )
  col_order <- col_order[names(col_order) %in% colnames(res)]
  setorderv(res, cols = names(col_order), order = col_order)

  setDF(res)

  return(res)
})


# Save
usethis::use_data(INCIDENCE_DA, overwrite = TRUE, version = 3)
