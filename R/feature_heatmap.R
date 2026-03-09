feature_heatmap = function(interested_genes,
                           da_results = "",
                           scale = 1,
                           padj_cutoff = 0.05,
                           filter_nonsig = FALSE,
                           title = NULL) {
  da_results = da_results %>%
    filter(gene_symbol %in% interested_genes) %>%
    mutate(
      sex   = sub("^(.)_.*", "\\1", contrast),
      group = sub("^._(\\S+).*", "\\1", contrast)
    ) %>%
    group_by(contrast, gene_symbol) %>%
    slice_min(adj.P.Val) %>%
    ungroup()

  expected_groups = c("1W", "2W", "4W", "8W")
  da_res_tidy = da_results %>%
    group_by(sex) %>%
    tidyr::complete(gene_symbol = interested_genes,
                    group = expected_groups,
                    fill = list(z = NaN)) %>%
    ungroup() %>%
    select(gene_symbol, group, z, adj.P.Val, sex)

  make_matrix = function(value_col, group_levels = expected_groups) {
    da_res_tidy %>%
      mutate(group = factor(group, levels = group_levels)) %>%
      arrange(sex, group) %>%
      distinct(sex, group, gene_symbol, .keep_all = TRUE) %>%
      rename(value = !!value_col) %>%
      pivot_wider(
        id_cols = c(sex, group),
        values_from = value,
        names_from = gene_symbol
      ) %>%
      as.data.frame() %>%
      `rownames<-`(with(., interaction(sex, group))) %>%
      select(-c(sex, group)) %>%
      t() %>%
      as.matrix()
  }

  ome_metadata_z   = make_matrix("z")
  ome_metadata_fdr = make_matrix("adj.P.Val")

  if (filter_nonsig && nrow(ome_metadata_z) > 25) {
    sig_rows = apply(ome_metadata_fdr, 1, function(r) any(!is.na(r) & r < padj_cutoff))
    ome_metadata_z   = ome_metadata_z[sig_rows, , drop = FALSE]
    ome_metadata_fdr = ome_metadata_fdr[sig_rows, , drop = FALSE]
  }

  # Build annotation from column names (format: "F.1W", "M.1W", etc.)
  annotation_df = strsplit(colnames(ome_metadata_z), split = "\\.") %>%
    do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("sex", "group")) %>%
    mutate(
      Sex = ifelse(sex == "F", "Female", "Male"),
      Timepoint = factor(group, levels = expected_groups)
    )

  top_annotation = ComplexHeatmap::HeatmapAnnotation(
    df = annotation_df %>% select(Sex, Timepoint),
    border = TRUE,
    gp = grid::gpar(col = "black"),
    gap = grid::unit(0, "pt"),
    which = "column",
    height = grid::unit(6 * 2, "pt") * scale,
    col = list(
      Sex = c(Female = MotrpacBicQC::sex_cols[["Female"]],
              Male   = MotrpacBicQC::sex_cols[["Male"]]),
      Timepoint = c("1W" = "#F7FCB9", "2W" = "#ADDD8E",
                    "4W" = "#238443", "8W" = "#002612")
    ),
    annotation_name_gp = grid::gpar(fontsize = 7 * scale),
    annotation_legend_param = list(
      border = "black",
      labels_gp = grid::gpar(fontsize = 6.5 * scale),
      title_gp  = grid::gpar(fontsize = 7 * scale, fontface = "bold")
    )
  )

  z_range = range(ome_metadata_z, na.rm = TRUE)
  max_abs  = max(abs(z_range))

  ht = ComplexHeatmap::Heatmap(
    matrix = ome_metadata_z,
    col = circlize::colorRamp2(
      breaks = c(-max_abs, 0, max_abs),
      colors = c("#3366ff", "white", "darkred")
    ),
    cluster_columns = FALSE,
    cluster_rows    = FALSE,
    show_column_names = FALSE,
    top_annotation  = top_annotation,
    border = "black",
    row_names_gp = grid::gpar(fontsize = 7 * scale),
    height = nrow(ome_metadata_z) * grid::unit(5.5, "pt") * scale,
    width  = ncol(ome_metadata_z) * grid::unit(5.5, "pt") * scale,
    column_split = annotation_df$Sex,
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 5 * scale, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Z-Score",
      at = c(-round(max_abs), 0, round(max_abs)),
      title_gp  = grid::gpar(fontsize = 7 * scale, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 6 * scale),
      legend_height = 5 * scale * grid::unit(8, "pt"),
      border = "black"
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.rect(x = x, y = y, width = width, height = height,
                      gp = grid::gpar(col = "#555555", fill = NA))
      if (!is.na(ome_metadata_fdr[i, j]) && ome_metadata_fdr[i, j] < padj_cutoff) {
        gb = grid::textGrob("*")
        gb_w = grid::convertWidth(grid::grobWidth(gb), "mm")
        gb_h = grid::convertHeight(grid::grobHeight(gb), "mm")
        grid::grid.text("*", x, y - gb_h * 0.5 + gb_w * 0.4)
      }
    }
  )

  ht
}
