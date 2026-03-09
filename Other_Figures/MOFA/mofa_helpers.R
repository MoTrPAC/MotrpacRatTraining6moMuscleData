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

#' Build a binary (0/1) matrix from a named list of feature sets
#' Rows = sets, columns = unique features across all sets
build_binary_matrix = function(feature_sets) {
  all_features = feature_sets %>% unlist() %>% unique()
  mat = matrix(0L, nrow = length(feature_sets), ncol = length(all_features))
  rownames(mat) = names(feature_sets)
  colnames(mat) = all_features
  for (i in seq_along(feature_sets)) {
    mat[i, feature_sets[[i]]] = 1L
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
                               csv_path = NULL) {
  feature_sets = lapply(feature_sets, toupper)
  binary_matrix = build_binary_matrix(feature_sets)

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

  # Export enrichment p-values to CSV
  if (!is.null(csv_path)) {
    enrichment_to_csv(enrichment_up, enrichment_down, view, csv_path)
  }

  invisible(list(up = enrichment_up, down = enrichment_down))
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
