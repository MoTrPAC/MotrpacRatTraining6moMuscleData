plot_rn7_data = function (assay = NULL, tissue = NULL, feature_ID = NULL, feature = NULL,
                          title = NULL, add_gene_symbol = FALSE, facet_by_sex = FALSE,
                          scale_x_by_time = TRUE, return_data = FALSE, exclude_outliers = TRUE,
                          add_adj_p = FALSE, ...)
{
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
  FEATURE = curr_feature
  redundant_feature = FEATURE
  differential = TRUE
  training_reg = data.table::as.data.table(MotrpacRatTraining6moData::TRAINING_REGULATED_FEATURES)
  keep_looking = TRUE
  while (keep_looking) {
    if (FEATURE %in% training_reg[, feature])
      break
    differential = FALSE
    if (!ASSAY %in% c("METAB", "IMMUNO"))
      break
    if (FEATURE %in% MotrpacRatTraining6moData::REPEATED_FEATURES$feature) {
      FEATURE = MotrpacRatTraining6moData::REPEATED_FEATURES$new_feature[MotrpacRatTraining6moData::REPEATED_FEATURES$feature ==
                                                                           FEATURE]
      differential = TRUE
      FEATURE_ID = MotrpacRatTraining6moData::REPEATED_FEATURES$feature_ID[MotrpacRatTraining6moData::REPEATED_FEATURES$new_feature ==
                                                                             FEATURE[1]]
      break
    }
    if (!ASSAY == "METAB")
      break
    new_feature_id = unique(MotrpacRatTraining6moData::METAB_FEATURE_ID_MAP$feature_ID_metareg[MotrpacRatTraining6moData::METAB_FEATURE_ID_MAP$metabolite_name ==
                                                                                                 FEATURE_ID])
    new_feature_id = unique(new_feature_id[new_feature_id %in%
                                             training_reg[assay == ASSAY & tissue == TISSUE,
                                                          feature_ID]])
    if (length(new_feature_id) == 0)
      break
    FEATURE = unique(training_reg[feature_ID == new_feature_id &
                                    tissue == TISSUE & assay == ASSAY, feature])
    if (length(FEATURE) == 0)
      break
    differential = TRUE
    keep_looking = FALSE
  }
  if (differential) {
    if (exclude_outliers) {
      all_sample_level_data = data.table::as.data.table(MotrpacRatTraining6moData::TRAINING_REGULATED_NORM_DATA_NO_OUTLIERS)
    }
    else {
      all_sample_level_data = data.table::as.data.table(MotrpacRatTraining6moData::TRAINING_REGULATED_NORM_DATA)
    }
    sample_level_data = all_sample_level_data[feature %in%
                                                FEATURE]
  }
  else {
    message(sprintf("'%s' is not a training-regulated feature. Looking in all sample-level data.",
                    FEATURE))
    all_sample_level_data = data.table::as.data.table(load_sample_data(TISSUE,
                                                                       ASSAY, exclude_outliers = exclude_outliers, ...))
    if (nrow(all_sample_level_data) == 0) {
      warning(sprintf("Sample-level data for %s %s not found.",
                      ASSAY, TISSUE))
      return()
    }
    if (!FEATURE_ID %in% all_sample_level_data[, feature_ID]) {
      warning(sprintf("'%s' not found in the %s %s sample-level data.",
                      FEATURE_ID, ASSAY, TISSUE))
      return()
    }
    sample_level_data = all_sample_level_data[feature_ID ==
                                                FEATURE_ID]
    sample_level_data[is.na(feature), `:=`(feature, FEATURE)]
  }
  multiple_measurements = FALSE
  if (nrow(sample_level_data) > 1) {
    message(sprintf("Multiple features correspond to '%s'. Plotting them together.",
                    redundant_feature))
    sample_level_data[, `:=`(feature, dataset)]
    multiple_measurements = TRUE
  }
  if (add_gene_symbol) {
    if (ASSAY %in% c("METHYL", "ATAC") & !differential) {
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
      title = sprintf("%s (%s)", redundant_feature, gene_symbol)
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
  value_cols = colnames(sample_level_data)[grepl("^[0-9]",
                                                 colnames(sample_level_data))]
  melted_subset = data.table::melt(sample_level_data, id.vars = c("feature"),
                                   measure.vars = value_cols, variable.name = "sample")
  melted_subset = melted_subset[!is.na(value)]
  melted_subset[, `:=`(sample, as.character(sample))]
  meta = unique(data.table::as.data.table(MotrpacRatTraining6moData::PHENO[,
                                                                           c("group", "sex", "pid", "viallabel")]))
  meta[, `:=`(pid, as.character(pid))]
  if (all(melted_subset[, sample] %in% meta[, viallabel])) {
    col = "viallabel"
  }
  else if (all(melted_subset[, sample] %in% meta[, pid])) {
    col = "pid"
    meta[, `:=`(viallabel, NULL)]
    meta = unique(meta)
  }
  else {
    stop(sprintf("Sample names in sample-level data do not correspond to vial labels or PIDs: %s...",
                 paste(utils::head(melted_subset[, sample]), collapse = ", ")))
  }
  subset_meta = merge(melted_subset, meta, by.x = "sample",
                      by.y = col)
  bygroup = subset_meta[, list(expr = mean(value, na.rm = T),
                               sd = sd(value, na.rm = T)), by = .(sex, group, feature)]
  if (return_data) {
    return(as.data.frame(bygroup))
  }
  bygroup[, `:=`(plot_group, sprintf("%s_%s", feature, sex))]
  if (multiple_measurements) {
    if (!facet_by_sex) {
      g = ggplot2::ggplot(bygroup, ggplot2::aes(x = group,
                                                y = expr, group = plot_group, colour = sex)) +
        ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::geom_point(size = 2, aes(shape = feature),
                            position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = expr -
                                              sd, ymax = expr + sd), width = 0.2, position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::theme_classic() + ggplot2::scale_colour_manual(values = c(female = MotrpacRatTraining6moData::SEX_COLORS[["F"]],
                                                                           male = MotrpacRatTraining6moData::SEX_COLORS[["M"]])) +
        ggplot2::labs(x = "Time trained (weeks)", y = "Normalized value",
                      title = title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                                                         size = 11), legend.title = ggplot2::element_blank(),
                                                      plot.subtitle = ggplot2::element_text(hjust = 0.5),
                                                      legend.position = "top", panel.grid.major = ggplot2::element_blank(),
                                                      panel.grid.minor = ggplot2::element_blank(),
                                                      legend.margin = ggplot2::margin(t = -5, b = -5,
                                                                                      unit = "pt"), legend.spacing.y = ggplot2::unit(0,
                                                                                                                                     "pt"))
    }
    else {
      g = ggplot2::ggplot(bygroup, ggplot2::aes(x = group,
                                                y = expr, group = plot_group)) + ggplot2::geom_line(colour = MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]],
                                                                                                    position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::geom_point(colour = MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]],
                            aes(shape = feature), position = ggplot2::position_dodge(width = 0.3),
                            size = 2) + ggplot2::geom_errorbar(aes(ymin = expr -
                                                                     sd, ymax = expr + sd), width = 0.2, colour = MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]],
                                                               position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::theme_classic() + ggplot2::labs(x = "Time trained (weeks)",
                                                 y = "Normalized value", title = title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                                                                                                            size = 11), legend.title = ggplot2::element_blank(),
                                                                                                         plot.subtitle = ggplot2::element_text(hjust = 0.5),
                                                                                                         legend.position = "top", panel.grid.major = ggplot2::element_blank(),
                                                                                                         panel.grid.minor = ggplot2::element_blank(),
                                                                                                         legend.margin = ggplot2::margin(t = -5, b = -5,
                                                                                                                                         unit = "pt")) + ggplot2::facet_wrap(~sex)
    }
  }
  else {
    if (!facet_by_sex) {
      g = ggplot2::ggplot(bygroup, ggplot2::aes(x = group,
                                                y = expr, group = plot_group, colour = sex)) +
        ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = expr -
                                              sd, ymax = expr + sd), width = 0.2, position = ggplot2::position_dodge(width = 0.3)) +
        ggplot2::theme_classic() + ggplot2::scale_colour_manual(values = c(female = MotrpacRatTraining6moData::SEX_COLORS[["F"]],
                                                                           male = MotrpacRatTraining6moData::SEX_COLORS[["M"]])) +
        ggplot2::labs(x = "Time trained (weeks)", y = "Normalized value",
                      title = title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                                                         size = 11), legend.title = ggplot2::element_blank(),
                                                      plot.subtitle = ggplot2::element_text(hjust = 0.5),
                                                      legend.position = "top", panel.grid.major = ggplot2::element_blank(),
                                                      panel.grid.minor = ggplot2::element_blank(),
                                                      legend.margin = ggplot2::margin(t = -5, b = -5,
                                                                                      unit = "pt"), legend.spacing.y = ggplot2::unit(0,
                                                                                                                                     "pt"))
    }
    else {
      g = ggplot2::ggplot(bygroup, ggplot2::aes(x = group,
                                                y = expr, group = plot_group)) + ggplot2::geom_line(colour = MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]]) +
        ggplot2::geom_point(colour = MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]]) +
        ggplot2::geom_errorbar(aes(ymin = expr - sd,
                                   ymax = expr + sd), width = 0.2, colour = MotrpacRatTraining6moData::TISSUE_COLORS[[TISSUE]]) +
        ggplot2::theme_classic() + ggplot2::labs(x = "Time trained (weeks)",
                                                 y = "Normalized value", title = title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                                                                                                            size = 11), legend.title = ggplot2::element_blank(),
                                                                                                         plot.subtitle = ggplot2::element_text(hjust = 0.5),
                                                                                                         legend.position = "none", panel.grid.major = ggplot2::element_blank(),
                                                                                                         panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::facet_wrap(~sex)
    }
  }
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
  if (add_adj_p) {
    message("Adding differential analysis p-value...")
    if (use_feature) {
      da = plot_feature_logfc(feature = FEATURE, add_adj_p = TRUE,
                              return_data = TRUE)
    }
    else {
      da = plot_feature_logfc(assay = ASSAY, tissue = TISSUE,
                              feature_ID = original_feature_ID, add_adj_p = TRUE,
                              return_data = TRUE)
    }
    if (!is.null(da)) {
      adj_p_value = min(unique(da$selection_fdr), na.rm = TRUE)
      subtitle = sprintf("adj. p-value: %s", signif(adj_p_value,
                                                    digits = 2))
      g = g + labs(subtitle = subtitle)
    }
  }
  return(g)
}
