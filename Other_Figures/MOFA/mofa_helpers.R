# Helper functions for MOFA exploration
# Source this file at the top of Mofa_exploration.Rmd

#' Print a summary of an ExpressionSet (features, samples, design table)
describe_eset = function(eset) {
  cat("Features:", nrow(eset), "\n")
  cat("Samples:", ncol(eset), "\n")
  table(pData(eset)$sex, pData(eset)$timepoint)
}

#' Select top n features by row variance
select_top_var = function(mat, n) {
  vars = apply(mat, 1, var)
  top_idx = order(vars, decreasing = TRUE)[seq_len(n)]
  mat[top_idx, ]
}

#' Remap matrix columns from viallabel to pid
remap_to_pid = function(mat, eset) {
  pd = pData(eset)
  pid_map = setNames(as.character(pd$pid), pd$viallabel)
  colnames(mat) = pid_map[colnames(mat)]
  mat
}

#' Build binary matrix aligned to actual MOFA feature names for a view
#'
#' MOFA renames duplicated features across views by appending _viewname.
#' This builds the binary matrix using the real MOFA feature names as columns
#' but maps gene set membership via base names (suffix stripped).
#'
#' @param feature_sets Named list of gene sets (already toupper'd).
#' @param mofa_features Character vector of feature names from the MOFA object
#'   for a specific view (from features_names(mofa_trained)[[view]]).
#' @param view Character, view name used by MOFA for suffix stripping.
build_binary_matrix_for_view = function(feature_sets, mofa_features, view) {
  # Strip _view suffix that MOFA adds to duplicated features
  suffix = paste0("_", toupper(view), "$")
  base_names = sub(suffix, "", toupper(mofa_features))

  mat = matrix(0L, nrow = length(feature_sets), ncol = length(mofa_features))
  rownames(mat) = names(feature_sets)
  colnames(mat) = mofa_features
  for (i in seq_along(feature_sets)) {
    hits = mofa_features[base_names %in% feature_sets[[i]]]
    if (length(hits) > 0) mat[i, hits] = 1L
  }
  mat
}

#' Run MOFA enrichment (up and down), plot results per factor, and export CSV
#'
#' @param mofa_trained Trained MOFA object
#' @param view Character, name of the view
#' @param feature_sets Named list of feature sets (will be toupper'd)
#' @param factors Integer vector of factors to test
#' @param max_pathways Max pathways to show in each plot (default 15)
#' @param csv_path Path to write enrichment results CSV. If NULL, no CSV.
#'
#' @return Named list with "up" and "down" enrichment results (invisible)
run_mofa_enrichment = function(mofa_trained, view, feature_sets,
                               factors=1:10, max_pathways=15,
                               csv_path=NULL, plot=FALSE) {
  feature_sets = lapply(feature_sets, toupper)
  mofa_features = MOFA2::features_names(mofa_trained)[[view]]
  binary_matrix = build_binary_matrix_for_view(feature_sets, mofa_features, view)

  enrichment_down = run_enrichment(mofa_trained,
    view=view, factors=factors,
    feature.sets=binary_matrix,
    sign="negative",
    statistical.test="parametric"
  )

  enrichment_up = run_enrichment(mofa_trained,
    view=view, factors=factors,
    feature.sets=binary_matrix,
    sign="positive",
    statistical.test="parametric"
  )

  # Plot each factor, catching errors when no pathways are significant
  if (plot) {
    for (f in factors) {
      tryCatch(
        print(plot_enrichment(enrichment_up,
          factor=f, max.pathways=max_pathways, text_size=0.8
        )),
        error=function(e) cat("Factor", f, "(up): no significant pathways\n")
      )
      tryCatch(
        print(plot_enrichment(enrichment_down,
          factor=f, max.pathways=max_pathways, text_size=0.8
        )),
        error=function(e) cat("Factor", f, "(down): no significant pathways\n")
      )
    }
  }

  # Export enrichment p-values to CSV
  if (!is.null(csv_path)) {
    enrichment_to_csv(enrichment_up, enrichment_down, view, csv_path)
  }

  invisible(list(up = enrichment_up, down = enrichment_down))
}

#' Plot enrichment across all omes for a single factor (faceted by view)
#'
#' @param enrichment_list Named list of lists, each with "up" and "down" enrichment
#'   results from run_mofa_enrichment (names = view labels).
#' @param factor Integer, which factor to visualize.
#' @param max_pathways Max pathways per view (top by adjusted p-value).
#' @param padj_threshold Significance cutoff highlighted in the plot.
#' @param png_path Path for PDF export (despite the name, saves as PDF at 600
#'   dpi). Set to NA to skip export (default).
plot_enrichment_faceted = function(enrichment_list, factor, max_pathways = 8,
                                   padj_threshold = 0.05, png_path = NA,
                                   width = 16) {
  factor_col = paste0("Factor", factor)

  shorten_name = function(x, cap = 25) {
    ifelse(nchar(x) > cap, paste0(substr(x, 1, cap), "..."), x)
  }

  df_list = lapply(names(enrichment_list), function(view_name) {
    res = enrichment_list[[view_name]]

    extract_top = function(enrichment, sign_label) {
      padj = enrichment$pval.adj
      if (!factor_col %in% colnames(padj)) return(NULL)
      vals = padj[, factor_col]
      top_idx = order(vals)[seq_len(min(max_pathways, length(vals)))]
      data.frame(
        pathway = rownames(padj)[top_idx],
        padj = vals[top_idx],
        sign = sign_label,
        view = view_name,
        stringsAsFactors = FALSE
      )
    }

    rbind(
      extract_top(res$up, "positive"),
      extract_top(res$down, "negative")
    )
  })

  df = do.call(rbind, df_list)
  df$neg_log_padj = -log10(pmax(df$padj, 1e-10))
  df$neg_log_padj = ifelse(df$sign == "negative", -df$neg_log_padj, df$neg_log_padj)
  df$significant = df$padj < padj_threshold
  df$view = factor(df$view, levels=names(enrichment_list))
  # within each view, order pathways by score for readability
  df$pathway_label = shorten_name(df$pathway)
  # For duplicate shortened names, keep the most significant (lowest padj)
  df = df[order(df$padj), ]
  df = df[!duplicated(paste0(df$view, "|", df$pathway_label)), ]
  df$pathway_label = factor(df$pathway_label,
    levels=rev(unique(df$pathway_label[order(df$view, df$neg_log_padj)])))

  threshold_line = -log10(padj_threshold)
  p = ggplot2::ggplot(df, ggplot2::aes(x=neg_log_padj, y=pathway_label, fill=significant)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept=c(-threshold_line, threshold_line),
                        linetype="dashed", color="gray50") +
    ggplot2::facet_wrap(~view, scales="free_y", nrow=1) +
    ggplot2::scale_fill_manual(
      values=c("FALSE"="gray75", "TRUE"="#2166ac"),
      name=paste0("padj < ", padj_threshold)
    ) +
    ggplot2::labs(
      title=paste0("Factor ", factor, " enrichments"),
      x=expression(-log[10](padj) ~ "[pos = up-weighted, neg = down-weighted]"),
      y=NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text=ggplot2::element_text(face="bold"),
      axis.text.y=ggplot2::element_text(size=7)
    )

  if (!is.na(png_path)) {
    out_path = sub("\\.(png|pdf)$", ".png", png_path, ignore.case = TRUE)
    ggplot2::ggsave(out_path, plot = p, dpi = 600, width = width, height = 4, units = "in")
    cat("Saved:", out_path, "\n")
  }

  p
}

#' Load and combine significant enrichment results from all omes
#'
#' Reads all *_sig.csv files from the saved_outputs directory and returns a
#' single data.frame sorted by factor number then padj.
#'
#' @param dir Path to the saved_outputs directory.
#'
#' @return A data.frame with columns: feature_set, factor, pval, padj, sign,
#'   view, factor_num (integer for sorting).
load_sig_enrichments = function(dir = "saved_outputs") {
  sig_files = list.files(dir, pattern = "_sig\\.csv$", full.names = TRUE)
  if (length(sig_files) == 0) stop("No *_sig.csv files found in: ", dir)

  df = do.call(rbind, lapply(sig_files, read.csv, stringsAsFactors = FALSE))
  df$factor_num = as.integer(sub("Factor", "", df$factor))
  df = df[order(df$factor_num, df$view, df$padj), ]
  rownames(df) = NULL
  df
}

#' Summarize significant enrichments as a nested list by factor
#'
#' @param sig_df Data.frame from load_sig_enrichments().
#' @param padj_threshold Significance threshold (default 0.05).
#'
#' @return Named list keyed by factor (e.g. "Factor1"), each element a
#'   data.frame of significant pathways across all omes for that factor.
enrichments_by_factor = function(sig_df, padj_threshold = 0.05) {
  sig_df = sig_df[sig_df$padj < padj_threshold, ]
  factors = unique(sig_df$factor[order(sig_df$factor_num)])
  out = lapply(factors, function(f) {
    sub_df = sig_df[sig_df$factor == f, ]
    sub_df[order(sub_df$view, sub_df$sign, sub_df$padj), ]
  })
  names(out) = factors
  out
}

#' Plot individual feature weightings for selected factors, faceted by view
#'
#' Extracts the top \code{n_features} features (by absolute weight) from each
#' view for each requested factor and draws a lollipop/dot plot.  Positive
#' weights point right and negative weights point left so the sign is
#' immediately readable.
#'
#' @param mofa_trained Trained MOFA2 object.
#' @param factors Integer vector of factors to display (default 1:3).
#' @param views Character vector of view names to include.  NULL (default) uses
#'   all views.
#' @param n_features Number of top features per view per factor (by |weight|).
#' @param view_colors Named character vector mapping view names to colors.  NULL
#'   uses a built-in palette.
#' @param png_path Path for PNG export at 600 dpi.  NA (default) skips export.
#'
#' @return A ggplot2 object (invisibly when \code{png_path} is set).
plot_feature_weights = function(mofa_trained, factors = 1:3, views = NULL,
                                n_features = 10, view_colors = NULL,
                                png_path = NA) {
  all_views = MOFA2::views_names(mofa_trained)
  if (is.null(views)) views = all_views
  views = intersect(views, all_views)
  if (length(views) == 0) stop("No matching views found in the MOFA object.")

  weights_list = MOFA2::get_weights(mofa_trained, views = views, as.data.frame = TRUE)
  # get_weights returns a data.frame with columns: feature, factor, view, value
  factor_labels = paste0("Factor", factors)
  df = weights_list %>%
    filter(factor %in% factor_labels, view %in% views) %>%
    group_by(factor, view) %>%
    arrange(desc(abs(value))) %>%
    slice_head(n = n_features) %>%
    ungroup() %>%
    mutate(
      feature = sub(paste0("_(",  paste(toupper(views), collapse="|"), ")$"), "",
                    feature, ignore.case = TRUE),
      factor_num = as.integer(sub("Factor", "", factor)),
      factor_label = paste0("Factor ", factor_num),
      factor_label = factor(factor_label,
                            levels = paste0("Factor ", sort(unique(factor_num)))),
      view = factor(view, levels = views)
    )

  # Order features within each facet by weight for readability
  df = df %>%
    arrange(factor_label, view, value) %>%
    mutate(feature_order = paste0(factor_label, "__", view, "__", feature),
           feature_order = factor(feature_order, levels = unique(feature_order)))

  if (is.null(view_colors)) {
    default_pal = c(
      Transcriptomics   = "#377EB8",
      Proteomics        = "#4DAF4A",
      Phosphoproteomics = "#FF7F00",
      Redox             = "#74add1",
      Metabolomics      = "#E41A1C"
    )
    view_colors = default_pal[intersect(names(default_pal), views)]
    missing = setdiff(views, names(view_colors))
    if (length(missing) > 0) {
      extra = scales::hue_pal()(length(missing))
      names(extra) = missing
      view_colors = c(view_colors, extra)
    }
  }

  p = ggplot2::ggplot(df, ggplot2::aes(x = value, y = feature_order, color = view)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = value, y = feature_order, yend = feature_order),
      linewidth = 0.4
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "gray40",
                        linewidth = 0.3) +
    ggplot2::facet_grid(view ~ factor_label, scales = "free_y", space = "free_y") +
    ggplot2::scale_color_manual(values = view_colors, guide = "none") +
    ggplot2::scale_y_discrete(labels = function(x) {
      lbl <- sub(".*__", "", x)
      ifelse(nchar(lbl) > 18, paste0(substr(lbl, 1, 18), "..."), lbl)
    }) +
    ggplot2::labs(
      title = "Top feature weights per factor",
      x = "Weight",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(face = "bold", size = 18),
      strip.text.y = ggplot2::element_text(face = "bold", angle = 0, size = 18),
      axis.text.y  = ggplot2::element_text(size = 16),
      axis.text.x  = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(size = 17),
      plot.title   = ggplot2::element_text(size = 19, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank()
    )

  if (!is.na(png_path)) {
    n_cols = length(factors)
    n_rows = length(views)
    ggplot2::ggsave(png_path, plot = p, dpi = 600,
                    width = 10, height = 1.5 + 0.2 * n_features * n_rows,
                    units = "in")
    cat("Saved:", png_path, "\n")
    invisible(p)
  } else {
    p
  }
}

#' Convert MOFA enrichment results to a data.frame and write to CSV
enrichment_to_csv = function(enrichment_up, enrichment_down, view,
                             csv_path) {
  extract_df = function(enrichment, sign) {
    # enrichment$pval is a matrix: pathways x factors
    pval = enrichment$pval
    padj = enrichment$pval.adj
    df = expand.grid(
      feature_set = rownames(pval),
      factor = colnames(pval),
      stringsAsFactors = FALSE
    )
    df$pval = as.vector(pval)
    df$padj = as.vector(padj)
    df$sign = sign
    df$view = view
    df
  }

  combined = rbind(
    extract_df(enrichment_up, "positive"),
    extract_df(enrichment_down, "negative")
  )
  combined = combined[order(combined$factor, combined$padj), ]

  write.csv(combined, csv_path, row.names = FALSE)
  cat("Enrichment results written to", csv_path, "\n")

  sig_path = sub("\\.csv$", "_sig.csv", csv_path)
  sig = combined[combined$padj < 0.05, ]
  write.csv(sig, sig_path, row.names = FALSE)
  cat("Significant results (padj < 0.05):", nrow(sig), "rows ->", sig_path, "\n")
}

#' Bubble heatmap of MOFA enrichment results for one view
#'
#' Rows = top significant gene sets, columns = factors.  Circle size encodes
#' -log10(min(padj_up, padj_down)) (capped at 5 mm); circle color encodes the
#' signed -log10(padj) (positive = up-weighted, negative = down-weighted).
#' A grey background rect marks significant cells (padj < padj_threshold).
#' Left annotation shows the view name.  Mirrors the style of
#' \code{make_camera_heatmap} from cmeans_helpers.R.
#'
#' @param enrichment Named list with \code{up} and \code{down} elements, each
#'   an object returned by \code{MOFA2::run_enrichment}.
#' @param view_name Character label for the view (used in left annotation).
#' @param ome_cols Named character vector mapping view names to colours.
#' @param col_fun A \code{circlize::colorRamp2} function for the signed
#'   -log10(padj) colour scale (should be shared across all views).
#' @param n_top Max number of gene sets to show (default 8).
#' @param padj_threshold Significance threshold for row selection and grey rect
#'   background (default 0.05).
#' @param cell_size \code{grid::unit} controlling cell width and height.
#'
#' @return A \code{ComplexHeatmap::Heatmap} or NULL when no significant sets.
make_mofa_enrichment_heatmap = function(enrichment, view_name, ome_cols, col_fun,
                                        n_top=8, padj_threshold=0.05,
                                        cell_size=grid::unit(6, "mm")) {
  factor_cols = intersect(colnames(enrichment$up$pval.adj),
                          colnames(enrichment$down$pval.adj))
  padj_up = enrichment$up$pval.adj[, factor_cols, drop=FALSE]
  padj_down = enrichment$down$pval.adj[, factor_cols, drop=FALSE]

  # Select top n_top gene sets with any significant result, ordered by min padj
  min_padj_by_set = pmin(apply(padj_up, 1, min, na.rm=TRUE),
                         apply(padj_down, 1, min, na.rm=TRUE))
  sig_sets = names(sort(min_padj_by_set[min_padj_by_set < padj_threshold]))
  if (length(sig_sets) == 0) return(NULL)
  top_sets = sig_sets[seq_len(min(n_top, length(sig_sets)))]

  padj_up_top = padj_up[top_sets, , drop=FALSE]
  padj_down_top = padj_down[top_sets, , drop=FALSE]

  # Signed -log10(padj): positive = up-weighted, negative = down-weighted; capped at ±20
  signed_logp = ifelse(padj_up_top <= padj_down_top,
                       -log10(padj_up_top),
                       log10(padj_down_top))
  signed_logp = pmax(pmin(signed_logp, 20), -20)
  min_padj_top = pmin(padj_up_top, padj_down_top)

  # Shorten long gene set names
  rn = ifelse(nchar(rownames(signed_logp)) > 40,
              paste0(substr(rownames(signed_logp), 1, 40), "..."),
              rownames(signed_logp))
  rownames(signed_logp) = rn
  rownames(min_padj_top) = rn
  colnames(signed_logp) = sub("Factor", "F", colnames(signed_logp))
  colnames(min_padj_top) = colnames(signed_logp)

  ComplexHeatmap::Heatmap(
    matrix=matrix(NA_real_, nrow=nrow(signed_logp), ncol=ncol(signed_logp),
                  dimnames=dimnames(signed_logp)),
    width=cell_size * ncol(signed_logp),
    height=cell_size * nrow(signed_logp),
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    rect_gp=grid::gpar(fill=NA, col="grey85"),
    show_column_names=TRUE,
    column_names_rot=0,
    column_names_gp=grid::gpar(just="left", fontsize=9),
    column_names_side="top",
    show_row_names=TRUE,
    row_names_side="right",
    row_names_gp=grid::gpar(fontsize=8),
    show_heatmap_legend=FALSE,
    na_col="grey95",
    left_annotation=ComplexHeatmap::rowAnnotation(
      ome=ComplexHeatmap::anno_simple(
        rep(view_name, nrow(signed_logp)),
        col=ome_cols,
        border=TRUE
      ),
      show_annotation_name=FALSE,
      width=grid::unit(4, "mm")
    ),
    cell_fun=function(j, i, x, y, width, height, fill) {
      pv = min_padj_top[i, j]
      fc = signed_logp[i, j]
      if (!is.na(pv) && pv < padj_threshold)
        grid::grid.rect(x, y, width, height, gp=grid::gpar(fill="grey90", col=NA))
      if (!is.na(pv) && pv < 1) {
        pv_size = min(-log10(pv), 5)
        if (!is.na(fc))
          grid::grid.circle(x=x, y=y, r=grid::unit(pv_size, "mm") / 2,
                            gp=grid::gpar(fill=col_fun(fc), col=NA))
      }
    }
  )
}


get_mofa_da_df = function(mofa_trained,
                          factor_num,
                          view,
                          da_results,
                          n=15,
                          gene_map=NULL,
                          feat_col="featureName") {

  weights = MOFA2::get_weights(mofa_trained, views=view, factors = factor_num, as.data.frame=TRUE) %>%
    dplyr::arrange(dplyr::desc(abs(value))) %>%
    dplyr::slice_head(n=n) %>%
    dplyr::arrange(dplyr::desc(value)) %>%
    dplyr::select(mofa_feature=feature, weight=value)

  if (!is.null(gene_map)) {
    mapped = gene_map(weights$mofa_feature)  # named vector: name=mofa_feature, value=da_id
    feature_map = data.frame(
      mofa_feature = names(mapped),
      da_id = unname(mapped),
      stringsAsFactors = FALSE
    )
  } else {
    feature_map = data.frame(
      mofa_feature = weights$mofa_feature,
      da_id = weights$mofa_feature,
      stringsAsFactors = FALSE
    )
  }

  da_lookup = da_results %>%
    dplyr::rename(da_id=!!feat_col) %>%
    dplyr::mutate(da_id = toupper(as.character(da_id)))

  weights %>%
    dplyr::left_join(feature_map, by="mofa_feature") %>%
    dplyr::left_join(da_lookup, by="da_id") %>%
    select(mofa_feature, weight, contrast, z, adj.P.Val) %>%
    mutate(z = as.numeric(z),
           adj.P.Val = as.numeric(adj.P.Val)) %>%
    arrange(desc(weight))
}


#' Heatmap of top MOFA features using z-scored normalized input data
#'
#' For each view, extracts the top \code{n} features by absolute MOFA weight,
#' pulls the normalized input data stored in the MOFA object via
#' \code{get_data()}, z-scores each feature across all samples, then averages
#' z-scores within each sex x timepoint group. Columns are sex x timepoint
#' groups; rows are features colored by view.
#'
#' @param mofa_trained Trained MOFA2 object.
#' @param factor_num Integer factor index (e.g. 1 for Factor1).
#' @param views Character vector of view names to include. NULL (default) uses
#'   all views.
#' @param n Number of top features per view (default 15).
#' @param scale Scalar multiplier for font sizes and cell dimensions (default 1).
#'
#' @return A \code{ComplexHeatmap::Heatmap} object, or NULL (invisibly) when no
#'   view yields matching features.
mofa_feature_heatmap_zscore = function(mofa_trained,
                                factor_num,
                                views = NULL,
                                n = 15,
                                scale = 1) {
  factor_label = paste0("Factor", factor_num)
  all_views = MOFA2::views_names(mofa_trained)
  if (is.null(views)) views = all_views
  views = intersect(views, all_views)

  # View colours matching the palette in plot_feature_weights
  default_view_cols = c(
    Transcriptomics   = "#377EB8",
    Proteomics        = "#4DAF4A",
    Phosphoproteomics = "#FF7F00",
    Redox             = "#74add1",
    Metabolomics      = "#E41A1C"
  )

  # Group order: Female then Male, within sex by timepoint
  tp_levels = c("SED", "1W", "2W", "4W", "8W")
  meta = MOFA2::samples_metadata(mofa_trained) %>%
    dplyr::mutate(timepoint = factor(timepoint, levels = tp_levels))

  # One column per sex x timepoint group, in order
  group_order = meta %>%
    dplyr::distinct(sex, timepoint) %>%
    dplyr::arrange(sex, timepoint)

  # All weights for this factor across requested views
  weights_df = MOFA2::get_weights(mofa_trained, views = views,
                                  factors = factor_label, as.data.frame = TRUE)

  mofa_data = MOFA2::get_data(mofa_trained)
  expr_mats = list()
  view_labels = c()

  for (v in views) {
    top_features = weights_df %>%
      dplyr::filter(view == v) %>%
      dplyr::arrange(dplyr::desc(abs(value))) %>%
      dplyr::slice_head(n = n) %>%
      #so first pull by abs value -> then pull and show positive first.
      dplyr::arrange(dplyr::desc(value)) %>%
      dplyr::pull(feature)

    if (length(top_features) == 0) next

    view_mat = mofa_data[[v]][[1]]  # feature x sample matrix
    colnames(view_mat) = gsub("^group1\\.", "", colnames(view_mat))

    feat_present = intersect(top_features, rownames(view_mat))

    if (length(feat_present) == 0) {
      message("No matching features for ", v, " Factor ", factor_num, " - skipping")
      next
    }

    # Z-score each feature across all samples first
    mat = view_mat[feat_present, meta$sample, drop = FALSE]
    mat = t(scale(t(mat)))

    # Average z-scores within each sex x timepoint group
    avg_mat = do.call(cbind, lapply(seq_len(nrow(group_order)), function(i) {
      grp_samples = meta$sample[meta$sex == group_order$sex[i] &
                                  meta$timepoint == group_order$timepoint[i]]
      grp_samples = intersect(grp_samples, colnames(mat))
      rowMeans(mat[, grp_samples, drop = FALSE], na.rm = TRUE)
    }))
    colnames(avg_mat) = paste(group_order$sex, group_order$timepoint, sep = ".")

    # Strip MOFA's _VIEWNAME suffix added to duplicated features
    rownames(avg_mat) = sub(paste0("_", toupper(v), "$"), "", rownames(avg_mat))

    expr_mats[[v]] = avg_mat
    view_labels = c(view_labels, rep(v, nrow(avg_mat)))
  }

  if (length(expr_mats) == 0) return(invisible(NULL))

  combined = do.call(rbind, expr_mats)

  # Column annotation: one row per sex x timepoint group
  ann_df = group_order %>%
    dplyr::rename(Sex = sex, Timepoint = timepoint)

  top_ann = ComplexHeatmap::HeatmapAnnotation(
    df = ann_df,
    border = TRUE,
    gp = grid::gpar(col = "black"),
    gap = grid::unit(0, "pt"),
    which = "column",
    height = grid::unit(8 * 2, "pt") * scale,
    col = list(
      Sex = c(Female = MotrpacBicQC::sex_cols[["Female"]],
              Male   = MotrpacBicQC::sex_cols[["Male"]]),
      Timepoint = c("SED" = "white", "1W" = "#F7FCB9", "2W" = "#ADDD8E",
                    "4W" = "#238443", "8W" = "#002612")
    ),
    annotation_name_gp = grid::gpar(fontsize = 6 * scale),
    annotation_legend_param = list(
      border = "black",
      labels_gp = grid::gpar(fontsize = 6.5 * scale),
      title_gp = grid::gpar(fontsize = 7 * scale, fontface = "bold")
    )
  )

  # Left annotation: view name per row
  present_views = unique(view_labels)
  view_cols = default_view_cols[intersect(names(default_view_cols), present_views)]
  missing_views = setdiff(present_views, names(view_cols))
  if (length(missing_views) > 0) {
    extra = scales::hue_pal()(length(missing_views))
    names(extra) = missing_views
    view_cols = c(view_cols, extra)
  }

  left_ann = ComplexHeatmap::rowAnnotation(
    View = ComplexHeatmap::anno_simple(view_labels, col = view_cols, border = TRUE),
    show_annotation_name = FALSE,
    width = grid::unit(4, "mm"),
    annotation_legend_param = list(
      title = "View",
      border = "black",
      labels_gp = grid::gpar(fontsize = 5.5 * scale),
      title_gp = grid::gpar(fontsize = 7 * scale, fontface = "bold")
    )
  )

  max_abs = max(abs(combined), na.rm = TRUE)

  ComplexHeatmap::Heatmap(
    matrix = combined,
    col = circlize::colorRamp2(c(-max_abs, 0, max_abs), c("#3366ff", "white", "darkred")),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    top_annotation = top_ann,
    left_annotation = left_ann,
    border = "black",
    row_names_gp = grid::gpar(fontsize = 6 * scale),
    height = nrow(combined) * grid::unit(6.5, "pt") * scale,
    width = ncol(combined) * grid::unit(5.5, "pt") * scale,
    column_split = ann_df$Sex,
    row_split = factor(view_labels, levels = names(expr_mats)),
    na_col = "grey90",
    heatmap_legend_param = list(
      title = "Scaled Expression",
      at = c(-round(max_abs), 0, round(max_abs)),
      title_gp = grid::gpar(fontsize = 7 * scale, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 6 * scale),
      legend_height = 5 * scale * grid::unit(8, "pt"),
      border = "black"
    )
  )
}


#' Heatmap of top MOFA features using differential analysis z-scores
#'
#' For each view in \code{view_da}, extracts the top \code{n} features by
#' absolute MOFA weight, maps them to DA results, and builds a heatmap where
#' columns are sex x timepoint contrasts (trained vs. sedentary) and cell color
#' encodes the limma z-score. Asterisks mark features with adj.P.Val <
#' \code{padj_cutoff}.
#'
#' @param mofa_trained Trained MOFA2 object.
#' @param factor_num Integer factor index (e.g. 1 for Factor1).
#' @param view_da Named list of DA result data frames (one per view), each with
#'   columns \code{gene_symbol}, \code{contrast}, \code{z}, \code{adj.P.Val}.
#' @param view_gene_map Named list of optional mapping functions
#'   \code{function(features) -> feature_ids}, one per view. NULL entries use
#'   direct matching.
#' @param n Number of top features per view (default 15).
#' @param scale Scalar multiplier for font sizes and cell dimensions (default 1).
#' @param padj_cutoff Significance threshold for asterisk annotation (default 0.05).
#' @param view_feature_col Named list specifying which fData column to match
#'   against MOFA feature names per view (default "featureName").
#'
#' @return A \code{ComplexHeatmap::Heatmap} object, or NULL (invisibly) when no
#'   view yields matching features.
mofa_feature_heatmap_da = function(mofa_trained,
                                   factor_num,
                                   view_da,
                                   view_gene_map = NULL,
                                   n = 15,
                                   scale = 1,
                                   padj_cutoff = 0.05,
                                   view_feature_col = NULL) {
  factor_label = paste0("Factor", factor_num)
  expected_groups = c("1W", "2W", "4W", "8W")

  default_view_cols = c(
    Transcriptomics   = "#377EB8",
    Proteomics        = "#4DAF4A",
    Phosphoproteomics = "#FF7F00",
    Redox             = "#74add1",
    Metabolomics      = "#E41A1C"
  )

  z_mats = list()
  fdr_mats = list()
  view_labels = c()

  for (v in names(view_da)) {
    gene_map = if (!is.null(view_gene_map)) view_gene_map[[v]] else NULL
    feat_col = if (!is.null(view_feature_col[[v]])) view_feature_col[[v]] else "featureName"

    da_df = get_mofa_da_df(mofa_trained, factor_num, v, view_da[[v]],
                           n = n, gene_map = gene_map, feat_col = feat_col)

    if (all(is.na(da_df$z))) {
      message("No matching features for ", v, " Factor ", factor_num, " - skipping")
      next
    }

    make_mat = function(value_col) {
      da_df %>%
        dplyr::filter(!is.na(contrast)) %>%
        dplyr::mutate(
          sex   = sub("^(.)_.*", "\\1", contrast),
          group = sub("^._(\\S+).*", "\\1", contrast)
        ) %>%
        dplyr::distinct(mofa_feature, sex, group, .keep_all = TRUE) %>%
        droplevels() %>%
        tidyr::pivot_wider(id_cols = mofa_feature,
                           names_from = c(sex, group),
                           names_sep = ".",
                           values_from = !!value_col) %>%
        tibble::column_to_rownames("mofa_feature") %>%
        as.matrix()
    }

    z_mats[[v]]   = make_mat("z")
    fdr_mats[[v]] = make_mat("adj.P.Val")
    view_labels = c(view_labels, rep(v, nrow(z_mats[[v]])))
  }

  if (length(z_mats) == 0) return(invisible(NULL))

  z_combined   = do.call(rbind, z_mats)
  fdr_combined = do.call(rbind, fdr_mats)

  ann_df = strsplit(colnames(z_combined), "\\.") %>%
    do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("sex", "group")) %>%
    dplyr::mutate(
      Sex      = ifelse(sex == "F", "Female", "Male"),
      Timepoint = factor(group, levels = expected_groups)
    )

  top_ann = ComplexHeatmap::HeatmapAnnotation(
    df = ann_df %>% dplyr::select(Sex, Timepoint),
    border = TRUE,
    gp = grid::gpar(col = "black"),
    gap = grid::unit(0, "pt"),
    which = "column",
    height = grid::unit(8 * 2, "pt") * scale,
    col = list(
      Sex = c(Female = MotrpacBicQC::sex_cols[["Female"]],
              Male   = MotrpacBicQC::sex_cols[["Male"]]),
      Timepoint = c("1W" = "#F7FCB9", "2W" = "#ADDD8E", "4W" = "#238443", "8W" = "#002612")
    ),
    annotation_name_gp = grid::gpar(fontsize = 7 * scale),
    annotation_legend_param = list(
      border = "black",
      labels_gp = grid::gpar(fontsize = 6.5 * scale),
      title_gp  = grid::gpar(fontsize = 7 * scale, fontface = "bold")
    )
  )

  present_views = unique(view_labels)
  view_cols = default_view_cols[intersect(names(default_view_cols), present_views)]
  missing_views = setdiff(present_views, names(view_cols))
  if (length(missing_views) > 0) {
    extra = scales::hue_pal()(length(missing_views))
    names(extra) = missing_views
    view_cols = c(view_cols, extra)
  }

  left_ann = ComplexHeatmap::rowAnnotation(
    View = ComplexHeatmap::anno_simple(view_labels, col = view_cols, border = TRUE),
    show_annotation_name = FALSE,
    width = grid::unit(4, "mm"),
    annotation_legend_param = list(
      title = "View",
      border = "black",
      labels_gp = grid::gpar(fontsize = 5.5 * scale),
      title_gp  = grid::gpar(fontsize = 7 * scale, fontface = "bold")
    )
  )

  z_range = range(z_combined, na.rm = TRUE)
  max_abs  = max(abs(z_range))

  ComplexHeatmap::Heatmap(
    matrix = z_combined,
    col = circlize::colorRamp2(c(-max_abs, 0, max_abs), c("#3366ff", "white", "darkred")),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    top_annotation = top_ann,
    left_annotation = left_ann,
    border = "black",
    row_names_gp = grid::gpar(fontsize = 6 * scale),
    height = nrow(z_combined) * grid::unit(6.5, "pt") * scale,
    width  = ncol(z_combined) * grid::unit(5.5, "pt") * scale,
    column_split = ann_df$Sex,
    row_split = factor(view_labels, levels = names(z_mats)),
    heatmap_legend_param = list(
      title = "DA Z-Score Relative\nto Sex-Matched Control",
      at = c(-round(max_abs), 0, round(max_abs)),
      title_gp  = grid::gpar(fontsize = 7 * scale, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 6 * scale),
      legend_height = 5 * scale * grid::unit(8, "pt"),
      border = "black"
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.rect(x = x, y = y, width = width, height = height,
                      gp = grid::gpar(col = "#555555", fill = NA))
      if (!is.na(fdr_combined[i, j]) && fdr_combined[i, j] < padj_cutoff) {
        gb   = grid::textGrob("*")
        gb_w = grid::convertWidth(grid::grobWidth(gb), "mm")
        gb_h = grid::convertHeight(grid::grobHeight(gb), "mm")
        grid::grid.text("*", x, y - gb_h * 0.5 + gb_w * 0.4)
      }
    }
  )
}



# Gene-level: most significant result per gene per contrast
get_gene_df = function(da_list, genes) {
  bind_rows(da_list) %>%
    mutate(row_id = toupper(gene_symbol)) %>%
    filter(!is.na(row_id), row_id %in% genes) %>%
    group_by(row_id, contrast) %>%
    slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(row_id, contrast, z, adj.P.Val)
}

# PTM site-level: match by gene name + position number in ptm_id.
# Uses position number (e.g. "195" from "S195") to handle both uppercase
# ("S195") and semicolon-prefixed (";195") ptm_id formats.
get_ptm_df = function(da_list, fdata, sites) {
  parsed = do.call(rbind, lapply(sites, function(s) {
    data.frame(
      gene = sub("_.*", "", s),
      pos = gsub("[^0-9]", "", sub(".*_", "", s)),
      row_id = s,
      stringsAsFactors = FALSE
    )
  }))

  fdata_df = fdata %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature_id") %>%
    mutate(gene_upper = toupper(gene_symbol))

  fmap = do.call(rbind, lapply(seq_len(nrow(parsed)), function(i) {
    pat = paste0("(^|[^0-9])", parsed$pos[i], "($|[^0-9])")
    fdata_df %>%
      filter(gene_upper == parsed$gene[i], grepl(pat, ptm_id)) %>%
      mutate(row_id = parsed$row_id[i]) %>%
      select(feature_id, row_id)
  }))

  if (is.null(fmap) || nrow(fmap) == 0) return(NULL)

  bind_rows(da_list) %>%
    inner_join(fmap, by = c("featureName" = "feature_id")) %>%
    group_by(row_id, contrast) %>%
    slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(row_id, contrast, z, adj.P.Val)
}

# Metabolite-level: match by metabolite_refmet (case-insensitive)
get_metab_df = function(da_list, metabolites) {
  target = toupper(metabolites)
  bind_rows(da_list) %>%
    mutate(row_id = toupper(metabolite_refmet)) %>%
    filter(!is.na(row_id), row_id %in% target) %>%
    group_by(row_id, contrast) %>%
    slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(row_id, contrast, z, adj.P.Val)
}
