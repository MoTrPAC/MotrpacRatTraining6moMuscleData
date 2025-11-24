library(MotrpacRatTraining6moMuscleData)
library(dplyr)
library(ComplexHeatmap)

#pull in all datasets into environment from zip file package
setwd("./data")
rdata_files <- list.files(path = ".", pattern = "\\.rda$", full.names = TRUE)
lapply(rdata_files, load, .GlobalEnv)


x <- list("GN-TRNSCRPT" = TRNSCRPT_GN_DA,
          "GN-PROT" = PROT_GN_DA,
          "GN-PHOSPHO" = PHOSPHO_GN_NORM_DA,
          "GN-ACETYL" = ACETYL_GN_NORM_DA,
          "GN-REDOX" = REDOX_GN_NORM_DA,
          "GN-METAB" = METAB_GN_DA,
          "VL-TRNSCRPT" = TRNSCRPT_VL_DA,
          "VL-METAB" = METAB_VL_DA) %>%
  lapply(function(x) {
    bind_rows(x) %>%
      summarise(.by = contrast,
                n_signif = sum(adj.P.Val < 0.05),
                n_total = sum(!is.na(logFC)))
  }) %>%
  bind_rows(.id = "ome") %>%
  mutate(prop_signif = n_signif / n_total)

prop_signif <- x %>%
  tidyr::pivot_wider(id_cols = ome,
                     names_from = contrast,
                     values_from = prop_signif) %>%
  as.data.frame() %>%
  `rownames<-`(.$ome) %>%
  select(-ome) %>%
  as.matrix()

n_signif <- x %>%
  tidyr::pivot_wider(id_cols = ome,
                     names_from = contrast,
                     values_from = n_signif) %>%
  as.data.frame() %>%
  `rownames<-`(.$ome) %>%
  select(-ome) %>%
  as.matrix()

rgb_to_greyscale <- function(col) {
  weights <- c("R" = 0.299, "G" = 0.587, "B" = 0.114)

  ceiling(apply(col2rgb(col) * weights, 2, sum)) / 255
}

color_fun <- circlize::colorRamp2(
  colors = c("#fcfdbfaa", "#fb8861aa", "#b63679aa", "#701090FF",
             rep("#100570ff", 2)),
  # breaks = 2 ^ (c(seq(0, log2(500 * ceiling(c(0.1, 0.5) * 100) / 100),
  #                   length.out = 10)) / 500)
  breaks = c(seq(0, 0.1, length.out = 5),
             ceiling(max(prop_signif) / 0.01) / 0.01)
)

layer_fun <- function(j, i, x, y, w, h, f) {
  grid.text(
    label = pindex(n_signif, i, j),
    x = x, y = y,
    gp = gpar(col = ifelse(
      col_fun(pindex(prop_signif, i, j)) %>%
        rgb_to_greyscale() %>%
        round(),
      "black",
      "white"
    ),
    fontsize = unit(9, "pt"))
  )
}


cell_size <- unit(14, "pt")

heatmap_args <- list(
  matrix = prop_signif,
  col = color_fun,
  border = TRUE,
  rect_gp = gpar(col = "grey40", lwd = 0.3),
  column_split = factor(rep(c("Female", "Male", ""), c(4, 4, 2)),
                        levels = c("Female", "Male", "")),
  heatmap_legend_param = list(
    title = "% DE Features (5% FDR):\n",
    title_position = "leftcenter",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    legend_direction = "horizontal",
    legend_width = unit(135, "pt"),
    # at = seq(0, 0.5, 0.1)
    at = seq(0, 0.10, 0.02),
    labels = latex2exp::TeX(
      c(as.character(100 * seq(0, 0.10 - 0.02, 0.02)), "$\\geq$10%")
    ),
    labels_gp = gpar(fontsize = 9)
  ),
  column_labels = latex2exp::TeX(
    c(rep(paste0(2^(0:3), "W - SED"), times = 2),
      "M$_{SED}$ - F$_{SED}$", "M$_{8W}$ - F$_{8W}$")
  ),
  # row_names_side = "left",
  row_names_gp = gpar(fontsize = 11),
  column_names_gp = gpar(fontsize = 11),
  left_annotation = HeatmapAnnotation(foo = anno_text(
    distinct(x, ome, n_total) %>%
      slice_max(n_total, by = ome) %>%
      pull(n_total) %>%
      {scales::number_format(big.mark = ",")(.)} %>%
      factor(., levels = .),
    gp = gpar(fontsize = 10),
    just = "right", location = 1),
    name = "",
    which = "row"
  ),
  layer_fun = layer_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  height = nrow(prop_signif) * cell_size,
  width = ncol(prop_signif) * cell_size * 1.7
)

col_fun <- heatmap_args$col

ht <- do.call(what = Heatmap, args = heatmap_args)

png("plots/SKM_DA_summary2.png",
    height = 3.5, width = 5.2, units = "in", res = 500)
draw(ht,
     heatmap_legend_side = "bottom")
dev.off()

