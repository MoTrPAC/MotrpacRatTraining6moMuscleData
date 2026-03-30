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
                               factors = 1:10, max_pathways = 15,
                               csv_path = NULL, plot = FALSE) {
  feature_sets = lapply(feature_sets, toupper)
  mofa_features = MOFA2::features_names(mofa_trained)[[view]]
  binary_matrix = build_binary_matrix_for_view(feature_sets, mofa_features, view)

  enrichment_down = run_enrichment(mofa_trained,
    view = view, factors = factors,
    feature.sets = binary_matrix,
    sign = "negative",
    statistical.test = "parametric"
  )

  enrichment_up = run_enrichment(mofa_trained,
    view = view, factors = factors,
    feature.sets = binary_matrix,
    sign = "positive",
    statistical.test = "parametric"
  )

  # Plot each factor, catching errors when no pathways are significant
  if (plot) {
    for (f in factors) {
      tryCatch(
        print(plot_enrichment(enrichment_up,
          factor = f, max.pathways = max_pathways, text_size = 0.8
        )),
        error = function(e) cat("Factor", f, "(up): no significant pathways\n")
      )
      tryCatch(
        print(plot_enrichment(enrichment_down,
          factor = f, max.pathways = max_pathways, text_size = 0.8
        )),
        error = function(e) cat("Factor", f, "(down): no significant pathways\n")
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
                                   padj_threshold = 0.05, png_path = NA) {
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
  df$view = factor(df$view, levels = names(enrichment_list))
  # within each view, order pathways by score for readability
  df$pathway_label = shorten_name(df$pathway)
  df$pathway_label = factor(df$pathway_label,
    levels = rev(unique(df$pathway_label[order(df$view, df$neg_log_padj)])))

  threshold_line = -log10(padj_threshold)
  p = ggplot2::ggplot(df, ggplot2::aes(x = neg_log_padj, y = pathway_label,
                                        fill = significant)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = c(-threshold_line, threshold_line),
                        linetype = "dashed", color = "gray50") +
    ggplot2::facet_wrap(~ view, scales = "free_y", nrow = 1) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "gray75", "TRUE" = "#2166ac"),
      name = paste0("padj < ", padj_threshold)
    ) +
    ggplot2::labs(
      title = paste0("Factor ", factor, " enrichment across all omes"),
      x = expression(-log[10](padj) ~ "[pos = up-weighted, neg = down-weighted]"),
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 7)
    )

  if (!is.na(png_path)) {
    out_path = sub("\\.(png|pdf)$", ".png", png_path, ignore.case = TRUE)
    ggplot2::ggsave(out_path, plot = p, dpi = 600, width = 16, height = 4, units = "in")
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
      Transcriptomics   = "#4575b4",
      Proteomics        = "#d73027",
      Phosphoproteomics = "#f46d43",
      Redox             = "#74add1",
      Metabolomics      = "#1a9850"
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
    ggplot2::scale_y_discrete(labels = function(x) sub(".*__", "", x)) +
    ggplot2::labs(
      title = "Top feature weights per factor",
      x = "Weight",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(face = "bold", size = 14),
      strip.text.y = ggplot2::element_text(face = "bold", angle = 0, size = 14),
      axis.text.y  = ggplot2::element_text(size = 13),
      axis.text.x  = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(size = 13),
      plot.title   = ggplot2::element_text(size = 15, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank()
    )

  if (!is.na(png_path)) {
    n_cols = length(factors)
    n_rows = length(views)
    ggplot2::ggsave(png_path, plot = p, dpi = 600,
                    width = 3.5 * n_cols, height = 1.5 + 0.25 * n_features * n_rows,
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
