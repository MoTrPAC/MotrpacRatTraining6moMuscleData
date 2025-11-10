## Purpose: Use WGCNA to determine clusters of related features

source("//protoapps/UserData/Sanford/For_Tyler/20210413_MoTrPAC_Data/Tyler_work/GSEA_functions.R")

library(writexl)
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
# disableWGCNAThreads()
library(RColorBrewer)
library(org.Mm.eg.db)
library(MSnSet.utils)
library(ComplexHeatmap)
library(glue)

# output folder
fs::dir_create("../output/WGCNA_global")

# process_data
# Load processed protein-level MSnSet
(load("../output/RD0a_1.1-msnsets_processed.RData"))
m3 <- m_ls[["m_glob"]]
# dim(exprs(m3))

# fix col name
fData(m3)$gene_name <- fData(m3)$gene_symbol

# Remove rows with missing gene symbols
m3 <- m3[!is.na(fData(m3)$gene_name), ]
# dim(exprs(m3))

# filter to 50% completeness
m3 <- filter_by_occurrence(m3, 0.5)

# This transpose is used multiple times throughout the script
datExpr <- t(exprs(m3))

## Network Construction and Module Detection ----
### Adjacency Matrix
# Choosing the soft power \beta
# The similarity matrix is constructed using
# biweight midcorrelation.
if (FALSE) {
    powers <- 5:30
    sft <- pickSoftThreshold(data = datExpr,
                             powerVector = powers,
                             networkType = "signed",
                             corFnc = WGCNA::bicor,
                             corOptions = list(use = "pairwise.complete.obs"),
                             verbose = 0,
                             RsquaredCut = 0.90)
    
    png(paste0(out_path, "/SoftPower_plot.png"), width = 10, height = 5, units = "in", res = 100)
    
    par(mfrow = c(1,2))
    cex1 = 0.9
    plot(sft$fitIndices[,1],
         -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"))
    text(sft$fitIndices[,1],
         -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red")
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.90,col="red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1],
         sft$fitIndices[,5],
         xlab="Soft Threshold (power)",
         ylab="Mean Connectivity",
         type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1],
         sft$fitIndices[,5],
         labels=powers,
         cex=cex1,
         col="red")
    dev.off()
}

# softPower <- sft$powerEstimate 
softPower <- 14

fs::dir_create(glue("../output/WGCNA_global/softPower_{softPower}"))
out_path <- paste0("../output/WGCNA_global/softPower_", softPower)


# Adjacency matrix
adjacency.mat <- adjacency(datExpr = datExpr,
                           power = softPower,
                           type = "signed",
                           corFnc = WGCNA::bicor,
                           corOptions = list(use = "pairwise.complete.obs"))

# A soft power of `r softPower` transforms the signed network so that
# it is approximately scale-free (so that the scale-free topology model
# fit $R^2$ is at least 0.90).

### Topological Overlap Matrix and Dissimilarity

# We compute the topological overlap matrix and its corresponding
# dissimilarity matrix (dissTOM = 1 - TOM).

# Topological Overlap Matrix (TOM)
TOM <- TOMsimilarity(adjacency.mat, TOMType = "unsigned")

# Dissimilarity matrix
dissTOM <- (1 - TOM)

### Average Link Hierarchical Clustering on dissTOM
# We use average link hierarchical clustering on the dissimilarity matrix
# to group similar proteins.

# Gene tree
geneTree <- hclust(as.dist(dissTOM), method = "average")

plot(geneTree,
     xlab = "", sub = "",
     ylab = "Dissimilarity",
     main = "Protein clustering on TOM-based dissimilarity",
     labels = FALSE,
     hang = 0.04)

#
# We will set the minimum number of proteins allowed in each cluster to 20.
# Reasonably numbers are around 20 or 30.
#

#+
# Dynamic tree cut
minModuleSize <- 20 # minimum number of genes per cluster
dynamicMods <- cutreeDynamic(dendro = geneTree,
                             distM = dissTOM,
                             minClusterSize = minModuleSize)

# length(unique(dynamicMods)) 

# Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods, colorSeq = standardColors())

#+ color_dendrogram_1, fig.wide = TRUE, fig.height = 6
# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Protein dendrogram and module colors")


# Modules can be merged if they are too similar.
# The similarity is quantified as correlation between the eigengenes.
# Threshold for module merged visualized on the dendrogram.
#

# Calculate eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
# ncol(MEs)

# Calculate dissimilarity of module eigengenes
MEDiss <- (1 - cor(MEs))

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# We choose a height cut of 0.2, which corresponds to a correlation of 0.8
MEDissThres <- 0.2

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h = MEDissThres, col = "red")

# Merge modules with correlation >= 0.8
merge <- mergeCloseModules(datExpr,
                           dynamicColors,
                           cutHeight = MEDissThres,
                           verbose = 3)
length(unique(merge$colors))

# merge modules
moduleColors <- merge$colors
# moduleColors <- dynamicColors

# Create a data frame with protein_name, Entry, gene_name, moduleColors,
# and sample abundances as columns
df.final <- cbind(fData(m3), moduleColors, exprs(m3))
# dim(df.final) 

write_xlsx(select(df.final, protein_id:moduleColors), path = glue::glue("{out_path}/RD1a_1-WGCNA_table_softpwr_{softPower}.xlsx"))


## ORA ----
# Map from gene symbol to Entrez ID
gene_clusters <- split(df.final$gene_name, df.final$moduleColors) %>%
  lapply(intersect, mappedRkeys(org.Mm.egSYMBOL2EG)) %>%
  lapply(function(mod_i) {
    sapply(as.list(org.Mm.egSYMBOL2EG[mod_i]), `[`, 1)
  }) %>%
  .[names(.) != "grey"] # Remove grey module


save(gene_clusters, file = glue::glue("{out_path}/RD1a_2-ORA_gene_clusters_sp{softPower}.RData"), compress = TRUE)

x <- compareCluster(geneClusters = gene_clusters,
                    fun = "enrichGO",
                    ont = "BP",
                    OrgDb = "org.Mm.eg.db",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    readable = TRUE,
                    universe = unlist(gene_clusters))

ora_results <- x


x@compareClusterResult %>%
  filter(p.adjust < 0.05) %>%
  arrange(Cluster, p.adjust) %>%
  write_xlsx(path = glue::glue("{out_path}/RD1a_3-WGCNA_ORA_BP_results_sp{softPower}.xlsx"))

res <- x@compareClusterResult %>%
  filter(p.adjust < 0.05) %>%
  group_by(Cluster) %>%
  arrange(p.adjust)

top_terms <- res %>%
  group_by(Cluster) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
  pull(ID) %>%
  unique()

res <- res %>%
  filter(ID %in% top_terms) %>%
  mutate(GeneProp = sapply(GeneRatio, function(generatio_i) {
    eval(parse(text = generatio_i))}),
    BgProp = sapply(BgRatio, function(generatio_i) {
      eval(parse(text = generatio_i))}),
    SetProp = Count / BgProp / 8172) %>%
  format_description() %>%
  arrange(p.adjust) %>%
  mutate(Description = factor(Description,
                              levels = unique(Description)),
         Cluster = factor(Cluster, levels = unique(Cluster)),
         p.adjust = ifelse(p.adjust <= 1e-06, 1e-06, p.adjust),
         p.adjust.discrete = cut(p.adjust, breaks = c(0, 1e-5, 0.0001, 0.001,
                                                      0.01, 0.05, 1),
                                 ordered_result = TRUE, right = FALSE))


p1 <- ggplot(res) +
  geom_point(aes(x = Cluster, y = Description, color = p.adjust),
             size = 2.7) +
  scale_x_discrete(name = "Module") +
  scale_y_discrete(name = NULL, limits = rev) +
  scale_color_viridis_c(
    name = "Adjusted\np-value",
    trans = "log10",
    breaks = trans_breaks("log10", function(x) 10 ^ x, n = 5),
    labels = label_scientific()
  ) +
  ggtitle("Top Overrepresented Biological Processes") +
  theme_bw(base_size = 8) +
  theme(text = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   color = "black", size = 7))

b <- 1.1
ggsave(glue::glue("{out_path}/RD1a_4-WGCNA_ORA_dotplot_pvalue_1_sp{softPower}.png"), p1,
       height = 7.8*b, width = 6.5*b, dpi = 300)


### Correlate with phenotypes ----

# eigengene msnset
ME_final <- merge$newMEs[, colnames(merge$newMEs) != "MEgrey"]
ME_final <- ME_final[sampleNames(m3),]
m_me <- MSnSet(exprs = t(as.matrix(ME_final)), pData = pData(m3))

ptype_to_test <- setdiff(varLabels(m_me), c("pid", "bid", "viallabel", "tmt11_channel", "tmt_plex",
                                            "timepoint", "exp_group"))

# res out 
lm_res <- vector("list", length(ptype_to_test))
names(lm_res) <- ptype_to_test

for (i in ptype_to_test) {
   lm_res[[i]] <- limma_gen(m_me, paste0("~ ", i), i) %>% 
      select(logFC, P.Value) %>% 
      rownames_to_column()
}


lm_df <- lm_res %>% 
   enframe() %>% 
   unnest(value) 

cor_mat <- lm_df %>% 
   select(-P.Value) %>% 
   pivot_wider(names_from = name, values_from = logFC) %>% 
   column_to_rownames() %>% 
   as.matrix()

p_mat <- lm_df %>% 
   select(-logFC) %>% 
   pivot_wider(names_from = name, values_from = P.Value) %>% 
   column_to_rownames() %>% 
   as.matrix() %>% 
   p.adjust(method = "BH") %>% 
   matrix(dim(cor_mat)[[1]], dim(cor_mat)[[2]])
colnames(p_mat) <- colnames(cor_mat)
rownames(p_mat) <- rownames(cor_mat)


# plot(y = exprs(m_me)[featureNames(m_me) == 'MEpink', ], x = pData(m_me)[, "delta_fat_pct"], asp = 1)
lm(exprs(m_me)[featureNames(m_me) == 'MEpink', ] ~ pData(m_me)[, "delta_fat_pct"]) %>% summary()
# abline(-0.0265, 0.03115)





# Markers for different levels of significance
comb_signif <- gtools::stars.pval(p_mat)

anno_colors <- list(module = sub("^ME", "", rownames(cor_mat)))
names(anno_colors$module) <- rownames(cor_mat)

# Legend for significance markers
gb <- textGrob(
  paste(".   < 0.1",
        "*   < 0.05",
        "**  < 0.01",
        "*** < 0.001",
        sep = "\n"),
  x = 0, y = 0,
  gp = gpar(fontsize = 6),
  hjust = 0, vjust = 0
)

lgd <- Legend(title = "Adjusted\np-value", grob = gb,
              # title_position = "lefttop-rot",
              title_gp = gpar(fontsize = 7))

# Column annotation
mod_count <- table(moduleColors)
mods <- names(mod_count)

# Determine color intensity for labels
pt_color <- sapply(mods, function(mod_i) {
  c( "black", "white")[
    1 + (sum(col2rgb(mod_i) * c(299, 587, 114)) / 1000 < 123)
  ]
})

ra <- HeatmapAnnotation(
  "Module Size" = anno_simple(rownames(cor_mat),
                              # border = TRUE,
                              pch = as.character(mod_count[rownames(mod_count) != "grey"]),
                              pt_gp = gpar(col = pt_color),
                              col = anno_colors$module,
                              width = unit(0.2, "in"),
                              pt_size = unit(6, "points")),
  which = "row", annotation_name_gp = gpar(fontsize = 7),
  show_legend = FALSE
)


ht <- Heatmap(cor_mat,
              column_title = "WGCNA Eigengene-Phenotype Correlation",
              cluster_columns = FALSE,
              border = TRUE,
              row_title_gp = gpar(fontsize = 8),
              column_title_gp = gpar(fontsize = 8),
              # show_heatmap_legend = FALSE,
              row_dend_side = "right",
              row_names_side = "left",
              column_title_side = "top",
              column_dend_height = unit(0.2, "in"),
              row_dend_width = unit(0.2, "in"),
              col = circlize::colorRamp2(breaks = c(-max(abs(cor_mat)), 0, max(abs(cor_mat))),
                                         colors = c("blue", "white", "red")),
              column_names_gp = gpar(fontsize = 7),
              row_names_gp = gpar(fontsize = 7),
              heatmap_legend_param = list(
                title = "logFC",
                title_gp = gpar(fontsize = 7),
                title_position = "lefttop-rot",
                legend_height = unit(0.85, "in"),
                grid_width = unit(0.1, "in"),
                labels_gp = gpar(fontsize = 6)
              ),
              left_annotation = ra,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(comb_signif[i, j], x, y,
                          gp = gpar(fontsize = 7))
              })

png(glue::glue("{out_path}/RD1a_5-WGCNA_module_cor_heatmap_sp{softPower}_limma.png"), 
    width = 4, 
    height = 3.5,
    units = "in", res = 300)
draw(ht, heatmap_legend_list = lgd, merge_legends = TRUE)
dev.off()
















