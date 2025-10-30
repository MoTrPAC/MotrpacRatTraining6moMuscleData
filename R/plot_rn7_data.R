#' Barebones version of `plot_feature_normalized_data` except
#' with rn7 data
#'
#' @param assay
#' @param tissue
#' @param feature_ID
#' @param feature
#' @param title
#' @param add_gene_symbol
#' @param facet_by_sex
#' @param scale_x_by_time
#' @param return_data
#' @param exclude_outliers
#' @param add_adj_p
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_rn7_data = function (assay = NULL,
                          tissue = NULL,
                          feature_ID = NULL,
                          feature = NULL,
                          title = NULL,
                          add_gene_symbol = FALSE,
                          scale_x_by_time = TRUE,
                          exclude_outliers = TRUE,
                          ...)
{
  assay <- match.arg(
    arg = assay,
    choices = c("TRNSCRPT"), #c("PROT", "TRNSCRPT", "METAB", "ACETYL", "METAB"),
    several.ok = FALSE
  )

  tissue <- match.arg(
    arg = tissue,
    choices = c("GN", "VL"),
    several.ok = FALSE
  )

  use_feature = !is.null(feature)
  curr_feature = feature
  if (is.null(curr_feature) & any(is.null(c(assay, tissue,
                                            feature_ID)))) {
    stop("If 'feature' is not specified, 'assay', 'tissue', and 'feature_ID' must all be specified.")
  }
  if (is.null(curr_feature)) {
    original_feature_ID = feature_ID
    if (assay == "METAB") {
      m = data.table(MotrpacRatTraining6moData::METAB_FEATURE_ID_MAP)
      if (!feature_ID %in% m[, feature_ID_sample_data]) {
        if (feature_ID %in% m[, metabolite_refmet]) {
          feature_ID = unique(m[metabolite_refmet ==
                                  feature_ID, feature_ID_sample_data])[1]
        }
        else {
          stop(sprintf("Feature ID '%s' not found in METAB data. See METAB_FEATURE_ID_MAP for measured metabolites.",
                       feature_ID))
        }
      }
    }
    curr_feature = sprintf("%s;%s;%s", assay, tissue, feature_ID)
  }
  if (is.null(tissue)) {
    splits = unname(unlist(strsplit(curr_feature, ";")))
    assay = splits[1]
    tissue = splits[2]
    feature_ID = splits[3]
    original_feature_ID = feature_ID
  }
  FEATURE_ID = feature_ID
  ASSAY = assay
  TISSUE = tissue

  #----trnscrpt only for the time being------
  obj_name = paste0(ASSAY, "_", TISSUE)
  data = get(obj_name, envir = asNamespace("MotrpacRatTraining6moMuscleData"))
  counts = exprs(data)
  dge = DGEList(counts = counts)
  dge = calcNormFactors(dge, method = "TMM")
  log2cpm = cpm(dge, log = TRUE, prior.count = 1)
  sample_level_data = log2cpm %>%
    as.data.frame() %>%
    filter(rownames(.) %in% FEATURE_ID) %>%
    log2() %>%
    t() %>%
    as.data.frame() %>%
    setNames("feature") %>%
    tibble::rownames_to_column("viallabel")
  #----trnscrpt only for the time being------

  if (add_gene_symbol) {
    if (ASSAY %in% c("METHYL", "ATAC")) {
      feature_to_gene = data.table::data.table(MotrpacRatTraining6moData::FEATURE_TO_GENE)
    }
    else {
      feature_to_gene = data.table::data.table(MotrpacRatTraining6moData::FEATURE_TO_GENE_FILT)
    }
    gene_symbol = feature_to_gene[feature_ID == FEATURE_ID,
                                  gene_symbol][1]
  }
  if (is.null(title)) {
    if (add_gene_symbol) {
      title = sprintf("%s", gene_symbol)
    }
    else {
      title = redundant_feature
    }
  }
  else {
    if (add_gene_symbol) {
      title = sprintf("%s (%s)", title, gene_symbol)
    }
  }
  subset_meta = MotrpacRatTraining6moData::PHENO %>%
    select(group, sex, pid, viallabel) %>%
    right_join(sample_level_data, by = "viallabel")

  bygroup = subset_meta %>%
    group_by(group) %>%
    summarise(
      expr = mean(feature, na.rm = TRUE),
      sd = sd(feature, na.rm = TRUE),
      .groups = "drop"
    )

  g = ggplot2::ggplot(bygroup, ggplot2::aes(x = group, y = expr, group = 1)) +
    ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.3)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.3)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = expr - sd, ymax = expr + sd),
      width = 0.2,
      position = ggplot2::position_dodge(width = 0.3)
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Time trained (weeks)",
      y = "Normalized value",
      title = title
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
      legend.position = "none",
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (scale_x_by_time) {
    g = g + scale_x_discrete(limits = c("control", "1w",
                                        "2w", "fill", "4w", rep("fill", 3), "8w"), labels = c("0",
                                                                                              "1", "2", "4", "8"), breaks = c("control", "1w",
                                                                                                                              "2w", "4w", "8w"))
  }
  else {
    g = g + scale_x_discrete(limits = c("control", "1w",
                                        "2w", "4w", "8w"), labels = c("0", "1", "2", "4",
                                                                      "8"), breaks = c("control", "1w", "2w", "4w", "8w"))
  }
  return(g)
}
