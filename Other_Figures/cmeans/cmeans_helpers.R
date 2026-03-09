# Helper functions for fuzzy c-means clustering analysis

# Reorder clusters so those with similar trajectories are consecutive.
reorder_clusters <- function(fclust) {
  d <- as.dist(1 - cor(t(fclust$centers)))
  hc <- hclust(d)
  neworder <- hc$order

  fclust$centers <- fclust$centers[neworder, ]
  rownames(fclust$centers) <- seq_len(nrow(fclust$centers))

  fclust$size <- fclust$size[neworder]

  fclust$cluster[] <- match(fclust$cluster, neworder)
  fclust$membership <- fclust$membership[, neworder]
  colnames(fclust$membership) <- seq_len(ncol(fclust$membership))

  return(fclust)
}

# Pivot a long data frame to a wide numeric matrix.
# id_col becomes row names, group_col becomes column names, value_col fills cells.
pivot_to_matrix <- function(df, id_col, group_col, value_col) {
  df %>%
    select(all_of(c(id_col, group_col, value_col))) %>%
    tidyr::pivot_wider(
      names_from = all_of(group_col),
      values_from = all_of(value_col)
    ) %>%
    as.data.frame() %>%
    `rownames<-`(.[[id_col]]) %>%
    select(-all_of(id_col)) %>%
    as.matrix()
}

# Create a colorRamp2 color function spanning 0 to the max of statistic_mat.
make_col_fun <- function(statistic_mat) {
  circlize::colorRamp2(
    c(0, max(statistic_mat, na.rm = TRUE)),
    c("white", "#503080")
  )
}

# Create a -log10(FDR) color legend using a given col_fun.
make_color_legend <- function(col_fun) {
  ComplexHeatmap::Legend(
    title = "-log10(FDR)",
    col_fun = col_fun,
    at = c(0, 10, 20),
    labels = c("0", "10", "20"),
    legend_height = grid::unit(4, "cm")
  )
}

# Build a bubble heatmap for CAMERA-PR results.
# Rows = gene sets, columns = sex_cluster_pattern groups.
# Bubble size encodes -log10(adj_p_value) (capped at 5 mm), color encodes statistic.
make_camera_heatmap <- function(statistic_mat, padj_mat, col_fun,
                                cell_size, sex_cols, ome_cols,
                                ome_i, is_first_ome) {
  sex_vec <- ifelse(
    startsWith(colnames(statistic_mat), "F"), "Female",
    ifelse(startsWith(colnames(statistic_mat), "M"), "Male", NA)
  )

  ComplexHeatmap::Heatmap(
    matrix = matrix(
      NA_real_,
      nrow = nrow(statistic_mat),
      ncol = ncol(statistic_mat),
      dimnames = dimnames(statistic_mat)
    ),
    width = cell_size * ncol(statistic_mat),
    height = cell_size * nrow(statistic_mat),

    cluster_rows = FALSE,
    cluster_columns = FALSE,

    rect_gp = grid::gpar(fill = NA, col = "grey85"),

    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = grid::gpar(just = "right"),
    show_row_names = TRUE,
    show_heatmap_legend = FALSE,
    row_names_side = "right",
    column_names_side = "top",
    na_col = "grey95",

    top_annotation =
      if (is_first_ome)
        ComplexHeatmap::HeatmapAnnotation(
          sex = ComplexHeatmap::anno_simple(sex_vec, col = sex_cols, border = FALSE),
          show_annotation_name = FALSE,
          height = grid::unit(4, "mm")
        )
      else NULL,

    left_annotation = ComplexHeatmap::rowAnnotation(
      ome = ComplexHeatmap::anno_simple(
        rep(ome_i, nrow(statistic_mat)),
        col = ome_cols,
        border = TRUE
      ),
      show_annotation_name = FALSE,
      width = grid::unit(4, "mm")
    ),
    heatmap_legend_param = list(
      at = c(0, 20),
      labels = c(0, 20)
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      fc <- statistic_mat[i, j]
      pv <- padj_mat[i, j]
      if (!is.na(pv) && pv < 0.05)
        grid::grid.rect(x, y, width, height, gp = grid::gpar(fill = "grey90", col = NA))
      if (!is.na(pv)) pv <- -log10(pv)
      if (!is.na(pv) && pv > 5) pv <- 5 # cap bubble size; colors are uncapped
      pv <- grid::unit(pv, "mm")
      if (!is.na(fc) && !is.na(pv)) {
        grid::grid.circle(
          x = x, y = y,
          r = pv / 2,
          gp = grid::gpar(fill = col_fun(fc), col = NA)
        )
      }
    }
  )
}
