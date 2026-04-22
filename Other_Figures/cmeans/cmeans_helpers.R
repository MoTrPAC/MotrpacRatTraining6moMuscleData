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

# Build a bubble heatmap for enrichment results, saving to filename.
# x:           pre-filtered data frame with columns set_col, cluster_num, adj_p_value, p_value
# set_col:     column name for gene set labels
# n_top:       top N terms per cluster (selected by smallest p_value among adj_p_value < 0.05)
# cluster_rows: whether to cluster heatmap rows
# height_extra / width_extra: extra inches added to the auto-computed dimensions
plot_enrichmap_clusters <- function(x, set_col = "GeneSet", n_top = 10L,
                                    cluster_rows = TRUE,
                                    height_extra = 1, width_extra = 6,
                                    filename) {
  ome_i <- if ("ome" %in% colnames(x)) unique(x$ome) else NULL
  padj_threshold <- 0.05 + 0.05 * isTRUE(length(ome_i) == 1L && ome_i == "PHOSPHO_GN")

  top_terms <- x %>%
    filter(adj_p_value < padj_threshold) %>%
    slice_min(p_value, n = n_top, by = cluster_num) %>%
    pull(all_of(set_col)) %>%
    unique()

  x <- x %>%
    filter(.data[[set_col]] %in% top_terms) %>%
    mutate(logP = -log10(p_value))

  height <- grid::convertUnit(length(top_terms) * grid::unit(15, "pt"), "in")
  height <- max(as.numeric(height) + 0.2, 4) + height_extra

  row_label_width <- grid::convertUnit(
    ComplexHeatmap::max_text_width(x[[set_col]]), "in"
  )
  width <- grid::convertUnit(nlevels(factor(x$cluster_num)) * grid::unit(15, "pt"), "in")
  width <- as.numeric(width + row_label_width) + width_extra

  enrichmap(
    x = x,
    plot_sig_only = FALSE,
    n_top = Inf,
    statistic_column = "logP",
    set_column = set_col,
    padj_column = "adj_p_value",
    padj_legend_title = "BH Adjusted\nP-Value",
    contrast_column = "cluster_num",
    heatmap_args = list(
      na_col = "grey95",
      rect_gp = grid::gpar(fill = "white", col = "grey85"),
      cluster_rows = cluster_rows,
      cluster_columns = FALSE,
      heatmap_legend_param = list(
        title = latex2exp::TeX("$\\bf{$-$log_{10}(P$-$Value)}$"),
        at = c(0, 20),
        labels = c(0, 20)
      ),
      column_names_side = "top"
    ),
    color = c("white", "#503080"),
    filename = filename,
    height = height,
    width = width
  )
}

# Run ORA (Fisher's exact) for each sex × ome × cluster combination.
# Returns a data frame with columns: sex, ome, cluster_num, set, p_value, adj_p_value, ...
run_fcm_ora <- function(eset, eset_list, FCM, index_list, num_clusters) {
  lapply(c("Female", "Male"), function(sex_i) {
    cl <- FCM[[sex_i]]
    eset_i <- eset_list[[sex_i]]

    lapply(unique(fData(eset)[["ome"]]), function(ome_i) {
      keep <- fData(eset_i)[["ome"]] == ome_i
      cluster <- cl$cluster[keep]
      names(cluster) <- fData(eset_i)[["featureName2"]][keep]
      cluster <- cluster[order(cluster)]
      cluster_list <- split(names(cluster),
                            factor(cluster, levels = as.character(unique(cluster))))

      background <- fData(eset_list[[sex_i]]) %>%
        filter(ome == ome_i) %>%
        pull(featureName2)

      lapply(seq_len(num_clusters), function(cluster_i) {
        res <- muscle_ora(
          input = cluster_list[[cluster_i]],
          background = background,
          database = index_list[[ome_i]]
        )
        res$cluster_num <- cluster_i
        res$adj_p_value <- p.adjust(res$p_value, method = "BH")
        res
      }) %>%
        bind_rows() %>%
        mutate(ome = ome_i)
    }) %>%
      bind_rows() %>%
      mutate(sex = sex_i)
  }) %>%
    bind_rows()
}

# Run a generic gene set test for each sex × ome combination.
# gst_fn(membership_mat, gene_sets) receives:
#   membership_mat: features × clusters numeric matrix (rownames = featureName2)
#   gene_sets:      index_list[[ome_i]]
# It should return a data frame; `ome` and `sex` columns are added automatically.
run_fcm_gst <- function(eset, eset_list, FCM, index_list, gst_fn) {
  lapply(c("Female", "Male"), function(sex_i) {
    cl <- FCM[[sex_i]]
    eset_i <- eset_list[[sex_i]]

    lapply(unique(fData(eset)[["ome"]]), function(ome_i) {
      keep <- fData(eset_i)[["ome"]] == ome_i
      membership_mat <- cl$membership[keep, ]
      rownames(membership_mat) <- fData(eset_i)[["featureName2"]][keep]

      gst_fn(membership_mat, index_list[[ome_i]]) %>%
        mutate(ome = ome_i)
    }) %>%
      bind_rows() %>%
      mutate(sex = sex_i)
  }) %>%
    bind_rows()
}
