# This documnet includes some figures that were previously around for mofa. No longer being used for the manuscript.
# ## Factor plots
#
# Visualize sample loadings on the first four factors, colored by experimental
# variables. I'm not sure why the plot_factors function is not working well for SEDS.
#
# Just plot through ggplot
#
# ```{r, fig.width=8, fig.height=6}
# factor_df = get_factors(mofa_trained, factors = 1:4, as.data.frame = TRUE) %>%
#   left_join(samples_metadata(mofa_trained), by = "sample") %>%
#   mutate(
#     factor = factor(factor, levels = paste0("Factor", 1:4)),
#     timepoint = factor(timepoint, levels = c("SED", "1W", "2W", "4W", "8W"))
#   )
#
# p_factor_sex = ggplot(factor_df, aes(x = "", y = value, fill = sex, shape = timepoint)) +
#   geom_jitter(width = 0.4, height = 0, size = 2.5, color = "black", stroke = 0.4) +
#   facet_wrap(~factor, nrow = 1) +
#   scale_shape_manual(values = c("SED" = 25, "1W" = 21, "2W" = 24, "4W" = 22, "8W" = 23)) +
#   scale_fill_manual(
#     values = c(Female = MotrpacBicQC::sex_cols[["Female"]],
#                Male = MotrpacBicQC::sex_cols[["Male"]])
#   ) +
#   labs(x = NULL, y = "Factor value", fill = "sex", shape = "timepoint") +
#   theme_bw() +
#   theme(
#     axis.text  = element_text(size = 12),
#     axis.title = element_text(size = 13),
#     strip.text = element_text(size = 12),
#     legend.text  = element_text(size = 11),
#     legend.title = element_text(size = 12)
#   )
# ggsave(file.path(here(), "plots", "MOFA", "factor_plot_by_sex.png"),
#        plot = p_factor_sex, width = 8, height = 6, units = "in", dpi = 600)
# p_factor_sex
# ```
#
# ```{r}
# #custom plots for panel 2, 3
# p = plot_feature_weights(mofa_trained,
#                          factors = 2,
#                          views = "Redox",
#                          n_features = 20,
#                          png_path = file.path(here(), "plots", "MOFA", "custom_panels", paste0("factor", 2, "_redox_feature_weights.png")))
#
# p = plot_feature_weights(mofa_trained,
#                          factors = 3,
#                          views = c("Proteomics",
#                                    "Phosphoproteomics",
#                                    "Redox"),
#                          n_features = 10,
#                          png_path = file.path(here(), "plots", "MOFA", "custom_panels", paste0("factor", 3, "_prot_feature_weights.png")))
# ```
#
# ```{r}
# ht = mofa_feature_heatmap_zscore(mofa_trained,
#                                  views = c("Proteomics",
#                                            "Phosphoproteomics",
#                                            "Redox"),
#                                  factor_num=3, n=8, scale=3)
#
# if (!is.null(ht)) {
#   png_path = file.path(here(), "plots", "MOFA", "custom_panels", paste0("factor", 3, "_prot_heatmap.png"))
#   png(png_path, width=10, height=12, units="in", res=600)
#   ComplexHeatmap::draw(ht, merge_legend=TRUE,
#                        column_title=paste0("Factor ", 3),
#                        column_title_gp=grid::gpar(fontsize=14, fontface="bold"))
#   dev.off()
#   cat("Saved:", png_path, "\n")
# }
#
# ```
#
#
# ```{r}
# phospho_mofa_to_da = function(mofa_features) {
#   m = site_map %>%
#     dplyr::filter(gene_site %in% mofa_features)
#   setNames(m$feature_id, m$gene_site)
# }
#
# redox_mofa_to_da = function(mofa_features) {
#   m = redox_feature_map %>%
#     dplyr::filter(gene_site %in% mofa_features)
#   setNames(m$feature_id, m$gene_site)
# }
#
# # deep dive specifically on Factor 2 redox
# ht = mofa_feature_heatmap_da(mofa_trained,
#                              factor_num=2,
#                              view_da=list(Redox = REDOX_GN_NORM_DA[["trained_vs_SED"]]),
#                              view_gene_map=list(Redox = redox_mofa_to_da),
#                              n=20,
#                              scale=3)
# if (!is.null(ht)) {
#   png_path = file.path(here(), "plots", "MOFA", "custom_panels", paste0("factor", 2, "_deepdive_heatmap.png"))
#   png(png_path, width=10, height=10, units="in", res=600)
#   ComplexHeatmap::draw(ht, merge_legend=TRUE,
#                        column_title=paste0("Factor ", 2),
#                        column_title_gp=grid::gpar(fontsize=14, fontface="bold"))
#   dev.off()
#   cat("Saved:", png_path, "\n")
# }
# ```
#
#
#
# ```{r}
# for (f in 1:10) {
#   p = plot_enrichment_faceted(all_enrichments,
#                               factor = f,
#                               max_pathways = 8,
#                               png_path = file.path(here(), "plots", "MOFA", paste0("factor", f, ".png")))
#   print(p)
# }
#
# ```
#
# Custom factor 2,3 plots
# ```{r}
# p = plot_enrichment_faceted(list(Redox = enrichment_redox),
#                             factor = 2,
#                             max_pathways = 7,
#                             width = 6,
#                             png_path = file.path(here(), "plots", "MOFA", "custom_panels", paste0("factor", 2, "_redox_pathways.png")))
# print(p)
#
# p = plot_enrichment_faceted(list(Proteomics = enrichment_prot,
#                                  Phosphoproteomics = enrichment_phospho,
#                                  Redox = enrichment_redox),
#                             factor = 3,
#                             max_pathways = 7,
#                             width = 12,
#                             png_path = file.path(here(), "plots", "MOFA","custom_panels",  paste0("factor", 3, "_protein_pathways.png")))
# print(p)
# ```
#
# ## Combined enrichment bubble heatmap
#
# Bubble heatmap across all views and factors.  Each row is a top significant
# gene set; columns are factors.  Circle size = -log10(padj) (capped at 5 mm),
# circle color = signed -log10(padj) (red = up-weighted, blue = down-weighted),
# grey background = padj < 0.05.
#
# ```{r mofa_enrichment_bubble, message=FALSE}
# # Compute global signed -log10(padj) range for a shared colour scale
# compute_signed_logp = function(enrich, n_top=8, padj_threshold=0.05) {
#   factor_cols = intersect(colnames(enrich$up$pval.adj),
#                           colnames(enrich$down$pval.adj))
#   pu = enrich$up$pval.adj[, factor_cols, drop=FALSE]
#   pd = enrich$down$pval.adj[, factor_cols, drop=FALSE]
#   min_p = pmin(apply(pu, 1, min, na.rm=TRUE), apply(pd, 1, min, na.rm=TRUE))
#   sig = names(sort(min_p[min_p < padj_threshold]))[seq_len(min(n_top, sum(min_p < padj_threshold)))]
#   if (length(sig) == 0) return(numeric(0))
#   as.vector(ifelse(pu[sig, , drop=FALSE] <= pd[sig, , drop=FALSE],
#                    -log10(pu[sig, , drop=FALSE]),
#                    log10(pd[sig, , drop=FALSE])))
# }
#
# all_vals = unlist(lapply(all_enrichments, compute_signed_logp))
# max_abs_logp = min(max(abs(all_vals), na.rm=TRUE), 20)
#
# col_fun_enrich = circlize::colorRamp2(
#   c(-max_abs_logp, 0, max_abs_logp),
#   c("#3366ff", "white", "darkred")
# )
#
# ome_cols_mofa = c(
#   Transcriptomics="#377EB8",
#   Proteomics="#4DAF4A",
#   Phosphoproteomics="#FF7F00",
#   Redox="#74add1",
#   Metabolomics="#E41A1C"
# )
#
# color_lg = ComplexHeatmap::Legend(
#   title="-log10(padj)\n(+ up, \u2212 down)",
#   col_fun=col_fun_enrich,
#   at=c(-round(max_abs_logp), 0, round(max_abs_logp)),
#   legend_height=grid::unit(4, "cm")
# )
# ome_lg = ComplexHeatmap::Legend(
#   title="View",
#   at=names(ome_cols_mofa),
#   legend_gp=grid::gpar(fill=ome_cols_mofa)
# )
# sig_lg = ComplexHeatmap::Legend(
#   title="Significance",
#   at="padj < 0.05",
#   legend_gp=grid::gpar(fill="grey90", col=NA)
# )
#
# enrich_heatmaps = list()
# for (v in names(all_enrichments)) {
#   ht = make_mofa_enrichment_heatmap(
#     enrichment=all_enrichments[[v]],
#     view_name=v,
#     ome_cols=ome_cols_mofa,
#     col_fun=col_fun_enrich,
#     n_top=8
#   )
#   if (!is.null(ht)) enrich_heatmaps[[v]] = ht
# }
#
# if (length(enrich_heatmaps) > 0) {
#   combined_enrich = Reduce(`%v%`, enrich_heatmaps)
#   png_path = file.path(here(), "plots", "MOFA", "enrichment_bubble_heatmap.png")
#   png(png_path, width=12, height=10, units="in", res=600)
#   ComplexHeatmap::draw(
#     combined_enrich,
#     heatmap_legend_list=list(color_lg, ome_lg, sig_lg),
#     merge_legend=TRUE,
#     heatmap_legend_side="left"
#   )
#   dev.off()
#   cat("Saved:", png_path, "\n")
# }
# ```
#
