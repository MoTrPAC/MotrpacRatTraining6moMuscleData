library(limma)
library(edgeR)
library(biomaRt)
library(PLIER)
library(ggplot2)
library(RColorBrewer)
library(qvalue)
library(genefilter)
library(gridExtra)
library(dplyr)
library(data.table)
library(ternvis)
library(ggpubr)
library(cowplot)
library(plyr)
library(stringr)
library(UpSetR)
library(foreign)
library(MASS)
library(sfsmisc)
library(IHW)
library("org.Rn.eg.db")
library(MCL)
library(corrplot)
library(ggrepel)
library(tidyverse)
library(ggsci)
library(circlize)
library(grid)
library(GenomicRanges)
library(readxl)
library(Gviz)
library(MotrpacHumanPreSuspension)

library(ChIPseeker)
require(TxDb.Rnorvegicus.UCSC.rn7.refGene)
txdb <- TxDb.Rnorvegicus.UCSC.rn7.refGene

library(Biobase)
library(edgeR)
library(ggvenn)


setwd("C:/Users/gsmit/Documents/Mount Sinai/Sealfon Laboratory/MoTrPAC/PASS1B Raw Data/LncRNA Project")

# This is the start of Gastroc PASS1B Analysis - Using the RN7 data

# We need the RN7 RNAseq Data for Gastroc, and RN7 ATACseq Data for Gastroc

# Loading the RN7 ATACseq Data
gastro_atacda <- read.delim(file = "RN7 ATAC DA and Norm Data 032825/pass1b-06_epigen-atac-seq_dea_rn7_pass1b-06_t55-gastrocnemius_epigen-atac-seq_training-dea-limma-ftest_20251703.txt",header = TRUE,sep = "\t",row.names = 1)
gastro_atacda_time <- read.delim(file = "RN7 ATAC DA and Norm Data 032825/pass1b-06_epigen-atac-seq_dea_rn7_pass1b-06_t55-gastrocnemius_epigen-atac-seq_timewise-dea-deseq2_20250704.txt",header = TRUE,sep = "\t")

gastro_atacda_time[is.na(gastro_atacda_time$adj_p_value),"adj_p_value"] <- 1

gastro_atacnorm <- read.delim(file = "RN7 ATAC DA and Norm Data 032825/pass1b-06_epigen-atac-seq_normalized-data_rn7_motrpac_pass1b-06_epigen-atac-seq_t55-gastrocnemius_limma_norm_20251703.tsv",header = TRUE,sep = "\t")

# Loading the RN7 RNAseq Data
load("TRNSCRPT_GN.rda")
load("TRNSCRPT_GN_DA.rda")

# Loading Homer-identified motif locations within each ATACseq peak
allpeakmotifs <- read.table(file = "pass1b_rn7_gastro_allpeak_allmotifs.txt",header = T,sep = "\t")

#gastro_atac_fw1sig <- gastro_atacda_time[gastro_atacda_time$sex %in% "female" & gastro_atacda_time$comparison_group %in% "1w" & gastro_atacda_time$adj_p_value < 0.1,"feature_ID"]
#gastro_atac_fw2sig <- gastro_atacda_time[gastro_atacda_time$sex %in% "female" & gastro_atacda_time$comparison_group %in% "2w" & gastro_atacda_time$adj_p_value < 0.1,"feature_ID"]
#gastro_atac_fw4sig <- gastro_atacda_time[gastro_atacda_time$sex %in% "female" & gastro_atacda_time$comparison_group %in% "4w" & gastro_atacda_time$adj_p_value < 0.1,"feature_ID"]
#gastro_atac_fw8sig <- gastro_atacda_time[gastro_atacda_time$sex %in% "female" & gastro_atacda_time$comparison_group %in% "8w" & gastro_atacda_time$adj_p_value < 0.1,"feature_ID"]

#gastro_atac_mw1sig <- gastro_atacda_time[gastro_atacda_time$sex %in% "male" & gastro_atacda_time$comparison_group %in% "1w" & gastro_atacda_time$adj_p_value < 0.1,"feature_ID"]
#gastro_atac_mw2sig <- gastro_atacda_time[gastro_atacda_time$sex %in% "male" & gastro_atacda_time$comparison_group %in% "2w" & gastro_atacda_time$adj_p_value < 0.1,"feature_ID"]
#gastro_atac_mw4sig <- gastro_atacda_time[gastro_atacda_time$sex %in% "male" & gastro_atacda_time$comparison_group %in% "4w" & gastro_atacda_time$adj_p_value < 0.1,"feature_ID"]
#gastro_atac_mw8sig <- gastro_atacda_time[gastro_atacda_time$sex %in% "male" & gastro_atacda_time$comparison_group %in% "8w" & gastro_atacda_time$adj_p_value < 0.1,"feature_ID"]

#gastro_atac_allsig <- Reduce(union,list(gastro_atac_fw1sig,
#                                        gastro_atac_fw2sig,
#                                        gastro_atac_fw4sig,
#                                        gastro_atac_fw8sig,
#                                        gastro_atac_mw1sig,
#                                        gastro_atac_mw2sig,
#                                        gastro_atac_mw4sig,
#                                        gastro_atac_mw8sig))

gastro_atac_training_sig <- rownames(gastro_atacda[gastro_atacda$p_value_adj < 0.1,])

gastro_atac_allmeta <- data.frame(row.names = rownames(gastro_atacnorm),
                                  "Chr" = gsub("chr","",gsub(":.*","",rownames(gastro_atacnorm))),
                                  "Start" = as.numeric(gsub(".*:","",gsub("-.*","",rownames(gastro_atacnorm)))),
                                  "End" = as.numeric(gsub(".*-","",rownames(gastro_atacnorm))))
gastro_atac_allmeta$Mid <- (gastro_atac_allmeta$Start + gastro_atac_allmeta$End)/2
gastro_atac_all_peak <- GRanges(seqnames = paste("chr",gastro_atac_allmeta$Chr,sep = ""), ranges = IRanges(gastro_atac_allmeta$Start, gastro_atac_allmeta$End))
gastro_atac_all_peakAnno <- annotatePeak(gastro_atac_all_peak,tssRegion=c(-2000,1000),TxDb = txdb,annoDb="org.Rn.eg.db")
gastro_atac_all_annogenes <- unique(gastro_atac_all_peakAnno@anno$ENSEMBL)
gastro_atac_all_peakAnnodf <- data.frame(row.names = rownames(gastro_atac_allmeta),
                                         "Chr" = gastro_atac_allmeta$Chr,
                                         "Start" = gastro_atac_allmeta$Start,
                                         "End" = gastro_atac_allmeta$End,
                                         "Mid" = gastro_atac_allmeta$Mid,
                                         "annotation" = gastro_atac_all_peakAnno@anno$annotation,
                                         "geneChr" = gastro_atac_all_peakAnno@anno$geneChr,
                                         "geneStart" = gastro_atac_all_peakAnno@anno$geneStart,
                                         "geneEnd" = gastro_atac_all_peakAnno@anno$geneEnd,
                                         "geneStrand" = gastro_atac_all_peakAnno@anno$geneStrand,
                                         "distanceToTSS" = gastro_atac_all_peakAnno@anno$distanceToTSS,
                                         "geneEnsembl" = gastro_atac_all_peakAnno@anno$ENSEMBL,
                                         "geneSymbol" = gastro_atac_all_peakAnno@anno$SYMBOL)

gastro_atac_allmeta$Peak <- rownames(gastro_atac_allmeta)

pa = as.data.table(gastro_atac_all_peakAnno@anno)

cols=c('seqnames','start','end')
gastropeaksdf <- as.data.table(gastro_atac_allmeta)
# This seems to do nothing
rownames(gastropeaksdf) <- gastro_atac_allmeta$Peak

if(nrow(pa)==nrow(gastropeaksdf)){
  pa[,feature_ID := gastropeaksdf[,"Peak"]]
}else{
  cols=c('seqnames','start','end')
  pa[,(cols) := lapply(.SD, as.character), .SDcols=cols]
  cols=c('chr','start','end')
  gastropeaksdf[,(cols) := lapply(.SD, as.character), .SDcols=cols]
  pa = merge(pa, gastropeaksdf, by.x=c('seqnames','start','end'), by.y=c('Chr','Start','End'), all.y=T)
}

pa[,short_annotation := annotation]
pa[grepl('Exon', short_annotation), short_annotation := 'Exon']
pa[grepl('Intron', short_annotation), short_annotation := 'Intron']

cols=c('start','end','geneStart','geneEnd','geneStrand')
pa[,(cols) := lapply(.SD, as.numeric), .SDcols=cols]
pa[,dist_upstream := ifelse(end-geneStart <=0, end-geneStart, NA_real_)]
pa[,dist_downstream := ifelse(start-geneEnd >= 0, start-geneEnd, NA_real_)]
# if feature overlaps with gene body, say dist_to_gene is 0
pa[end >= geneStart & start <= geneEnd, dist_downstream := 0]
pa[end >= geneStart & start <= geneEnd, dist_upstream := 0]
pa[, relationship_to_gene := ifelse(is.na(dist_downstream), dist_upstream, dist_downstream)]
pa[,c('dist_upstream','dist_downstream') := NULL]

pa[relationship_to_gene == 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Overlaps Gene"]
# Downstream
pa[geneStrand == 1 & relationship_to_gene > 0 & relationship_to_gene < 5000, short_annotation := "Downstream (<5kb)"]
pa[geneStrand == 2 & relationship_to_gene < 0 & relationship_to_gene > -5000, short_annotation := "Downstream (<5kb)"]
# Upstream (promoter excluded)
pa[geneStrand == 1 & relationship_to_gene > -5000 & relationship_to_gene < 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Upstream (<5kb)"]
pa[geneStrand == 2 & relationship_to_gene < 5000 & relationship_to_gene > 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Upstream (<5kb)"]
pa[abs(relationship_to_gene) >= 5000, short_annotation := "Distal Intergenic"]

setnames(pa, c('short_annotation','annotation','seqnames','geneId'), c('custom_annotation','chipseeker_annotation','chrom','ensembl_gene'))

gastro_atac_all_peakAnnodf <- data.frame(row.names = rownames(gastro_atac_allmeta),
                                         "Chr" = gastro_atac_allmeta$Chr,
                                         "Start" = gastro_atac_allmeta$Start,
                                         "End" = gastro_atac_allmeta$End,
                                         "Mid" = gastro_atac_allmeta$Mid,
                                         "annotation" = gastro_atac_all_peakAnno@anno$annotation,
                                         "geneChr" = gastro_atac_all_peakAnno@anno$geneChr,
                                         "geneStart" = gastro_atac_all_peakAnno@anno$geneStart,
                                         "geneEnd" = gastro_atac_all_peakAnno@anno$geneEnd,
                                         "geneStrand" = gastro_atac_all_peakAnno@anno$geneStrand,
                                         "distanceToTSS" = gastro_atac_all_peakAnno@anno$distanceToTSS,
                                         "geneEnsembl" = gastro_atac_all_peakAnno@anno$ENSEMBL,
                                         "geneSymbol" = gastro_atac_all_peakAnno@anno$SYMBOL,
                                         "custom_annotation" = pa$custom_annotation)

gastro_atac_sig_peakAnnodf <- gastro_atac_all_peakAnnodf[gastro_atac_training_sig,]



gastro_atac_sigl2fc <- matrix(0L,nrow = length(gastro_atac_training_sig),ncol = 8)
rownames(gastro_atac_sigl2fc) <- gastro_atac_training_sig
colnames(gastro_atac_sigl2fc) <- c("F W1","F W2","F W4","F W8","M W1","M W2","M W4","M W8")

for(i in 1:dim(gastro_atac_sigl2fc)[1]){
  if(i%%20 == 0){
    print(paste("gastro",toString(i),sep = " "))
  }
  ourpeak <- gastro_atac_training_sig[i]
  ourda <- gastro_atacda_time[gastro_atacda_time$feature_ID %in% ourpeak,]
  gastro_atac_sigl2fc[i,"F W1"] <- ourda[ourda$sex %in% "female" & ourda$comparison_group %in% "1w","logFC"]
  gastro_atac_sigl2fc[i,"F W2"] <- ourda[ourda$sex %in% "female" & ourda$comparison_group %in% "2w","logFC"]
  gastro_atac_sigl2fc[i,"F W4"] <- ourda[ourda$sex %in% "female" & ourda$comparison_group %in% "4w","logFC"]
  gastro_atac_sigl2fc[i,"F W8"] <- ourda[ourda$sex %in% "female" & ourda$comparison_group %in% "8w","logFC"]
  gastro_atac_sigl2fc[i,"M W1"] <- ourda[ourda$sex %in% "male" & ourda$comparison_group %in% "1w","logFC"]
  gastro_atac_sigl2fc[i,"M W2"] <- ourda[ourda$sex %in% "male" & ourda$comparison_group %in% "2w","logFC"]
  gastro_atac_sigl2fc[i,"M W4"] <- ourda[ourda$sex %in% "male" & ourda$comparison_group %in% "4w","logFC"]
  gastro_atac_sigl2fc[i,"M W8"] <- ourda[ourda$sex %in% "male" & ourda$comparison_group %in% "8w","logFC"]
  
}


# I want to do two things. First, I want an update TF motif-peak connection list.
#peaktableforhomer <- gastro_atac_all_peakAnnodf[,c("Chr","Start","End","geneStrand")]
#peaktableforhomer$Chr <- paste("chr",peaktableforhomer$Chr,sep = "")
#write.table(peaktableforhomer,file = "pass1bgastropeaktableforhomer.narrowPeak",sep = "\t",row.names = T,col.names = F,quote = F)

# I will eventually run PLIER on the ATAC data using associated gene pathways and TF enrichments as prior knowledge.


# First, let's take a closer look at the sig ATAC data.
hc <- hclust(dist(gastro_atac_sigl2fc), "complete")
ourclusters <- cutree(hc, k = 5)

# Figure S2A
png(file = "FigureS2A_gastro_atacsig_l2fc_heatmap_110725.png",width = 6,height = 6,units = "in",res = 600)
pheatmap(gastro_atac_sigl2fc,cluster_cols = F,show_rownames = F,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0)
dev.off()

pdf(file = "FigureS2A_gastro_atacsig_l2fc_heatmap_110725.pdf",width = 6,height = 6)
pheatmap(gastro_atac_sigl2fc,cluster_cols = F,show_rownames = F,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0)
dev.off()

#png(file = "gastro_atacsig_cluster1_l2fc_heatmap.png",width = 6,height = 6,units = "in",res = 600)
#pheatmap(gastro_atac_sigl2fc[names(ourclusters[ourclusters == 1]),],cluster_cols = F,show_rownames = T,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0)
#dev.off()

#png(file = "gastro_atacsig_cluster2_l2fc_heatmap.png",width = 6,height = 6,units = "in",res = 600)
#pheatmap(gastro_atac_sigl2fc[names(ourclusters[ourclusters == 2]),],cluster_cols = F,show_rownames = T,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0)
#dev.off()

#png(file = "gastro_atacsig_cluster3_l2fc_heatmap.png",width = 6,height = 6,units = "in",res = 600)
#pheatmap(gastro_atac_sigl2fc[names(ourclusters[ourclusters == 3]),],cluster_cols = F,show_rownames = T,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0)
#dev.off()

#png(file = "gastro_atacsig_cluster4_l2fc_heatmap.png",width = 6,height = 6,units = "in",res = 600)
#pheatmap(gastro_atac_sigl2fc[names(ourclusters[ourclusters == 4]),],cluster_cols = F,show_rownames = T,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0)
#dev.off()

#png(file = "gastro_atacsig_cluster5_l2fc_heatmap.png",width = 6,height = 6,units = "in",res = 600)
#pheatmap(gastro_atac_sigl2fc[names(ourclusters[ourclusters == 5]),],cluster_cols = F,show_rownames = T,breaks = seq(-3,3,length.out = 101),color = colorpanel(101,"blue","white","red"),angle_col = 0)
#dev.off()

cluster1sigpeaks <- names(ourclusters[ourclusters == 2])
cluster2sigpeaks <- names(ourclusters[ourclusters == 1])
cluster3sigpeaks <- names(ourclusters[ourclusters == 3])
cluster4sigpeaks <- names(ourclusters[ourclusters == 4])
cluster5sigpeaks <- names(ourclusters[ourclusters == 5])

clustersmedtraj <- rbind(apply(gastro_atac_sigl2fc[cluster1sigpeaks,],2,median),
                     apply(gastro_atac_sigl2fc[cluster2sigpeaks,],2,median),
                     apply(gastro_atac_sigl2fc[cluster3sigpeaks,],2,median),
                     apply(gastro_atac_sigl2fc[cluster4sigpeaks,],2,median),
                     apply(gastro_atac_sigl2fc[cluster5sigpeaks,],2,median))
rownames(clustersmedtraj) <- c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5")

clustersmeantraj <- rbind(apply(gastro_atac_sigl2fc[cluster1sigpeaks,],2,mean),
                         apply(gastro_atac_sigl2fc[cluster2sigpeaks,],2,mean),
                         apply(gastro_atac_sigl2fc[cluster3sigpeaks,],2,mean),
                         apply(gastro_atac_sigl2fc[cluster4sigpeaks,],2,mean),
                         apply(gastro_atac_sigl2fc[cluster5sigpeaks,],2,mean))
rownames(clustersmeantraj) <- c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5")

clustersmedtrajdf <- data.frame("L2FC" = as.vector(clustersmedtraj),
                                "Time Point" = rep(c("W1","W2","W4","W8",
                                                     "W1","W2","W4","W8"),each = 5),
                                "Cluster" = rep(c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5"),n = 8),
                                "Sex" = rep(c("F","F","F","F","M","M","M","M"),each = 5))

clustersmeantrajdf <- data.frame("L2FC" = as.vector(clustersmeantraj),
                                "Time Point" = rep(c("W1","W2","W4","W8",
                                                     "W1","W2","W4","W8"),each = 5),
                                "Cluster" = rep(c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5"),n = 8),
                                "Sex" = rep(c("F","F","F","F","M","M","M","M"),each = 5))

# Figure 2A
png("Figure2A_Gastro_ATAC_Cluster_Median_Trajectories_lineplot_110725.png",width = 6,height = 3,units = "in",res = 600)
ggline(clustersmedtrajdf,x = "Time.Point", y = "L2FC",color = "Cluster",facet.by = "Sex") + geom_hline(yintercept = 0,linetype = "dashed",colour = "grey75")
dev.off()

pdf("Figure2A_Gastro_ATAC_Cluster_Median_Trajectories_lineplot_110725.pdf",width = 6,height = 3)
ggline(clustersmedtrajdf,x = "Time.Point", y = "L2FC",color = "Cluster",facet.by = "Sex") + geom_hline(yintercept = 0,linetype = "dashed",colour = "grey75")
dev.off()

#png("Gastro_ATAC_Cluster_Mean_Trajectories_lineplot.png",width = 6,height = 3,units = "in",res = 600)
#ggline(clustersmeantrajdf,x = "Time.Point", y = "L2FC",color = "Cluster",facet.by = "Sex") + geom_hline(yintercept = 0,linetype = "dashed",colour = "grey75")
#dev.off()

# Let's make a table of annotation distributions for ATAC for all peaks, sig peaks and by cluster

gastro_atac_all_peakAnnodf$custom_annotation <- factor(gastro_atac_all_peakAnnodf$custom_annotation,levels = names(table(gastro_atac_all_peakAnnodf$custom_annotation)))

atacannotable <- cbind(table(gastro_atac_all_peakAnnodf$custom_annotation)/dim(gastro_atac_all_peakAnnodf)[1],
                       table(gastro_atac_all_peakAnnodf[gastro_atac_training_sig,"custom_annotation"])/length(gastro_atac_training_sig),
                       table(gastro_atac_all_peakAnnodf[cluster1sigpeaks,"custom_annotation"])/length(cluster1sigpeaks),
                       table(gastro_atac_all_peakAnnodf[cluster2sigpeaks,"custom_annotation"])/length(cluster2sigpeaks),
                       table(gastro_atac_all_peakAnnodf[cluster3sigpeaks,"custom_annotation"])/length(cluster3sigpeaks),
                       table(gastro_atac_all_peakAnnodf[cluster4sigpeaks,"custom_annotation"])/length(cluster4sigpeaks),
                       table(gastro_atac_all_peakAnnodf[cluster5sigpeaks,"custom_annotation"])/length(cluster5sigpeaks))
colnames(atacannotable) <- c("All Peaks","All Sig Peaks","Cluster1","Cluster2","Cluster3","Cluster4","Cluster5")

#Figure S2B
png("FigureS2B_Muscle ATAC Peak Distribution by Peak Type table.png",width = 6,height = 3,units = "in",res = 600)
pheatmap(atacannotable,cluster_rows = F,cluster_cols = F,angle_col = 315,color = colorpanel(101,"white","firebrick"),display_numbers = T,number_color = "black")
dev.off()

pdf("FigureS2B_Muscle ATAC Peak Distribution by Peak Type table.pdf",width = 6,height = 3)
pheatmap(atacannotable,cluster_rows = F,cluster_cols = F,angle_col = 315,color = colorpanel(101,"white","firebrick"),display_numbers = T,number_color = "black")
dev.off()

muscleatacsig_rat_to_human_df <- data.frame(row.names = gastro_atac_training_sig,
                                            "Peak" = gastro_atac_training_sig,
                                            "Rat_geneEnsembl" = gastro_atac_all_peakAnnodf[gastro_atac_training_sig,"geneEnsembl"])
muscleatacsig_rat_to_human_df$Human_geneEnsembl = ""
muscleatacsig_rat_to_human_df$Human_geneSymbol = ""
for(i in 1:length(gastro_atac_training_sig)){
  ourratens <- muscleatacsig_rat_to_human_df[i,"Rat_geneEnsembl"]
  if(ourratens %in% MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE$RAT_ENSEMBL_ID){
    muscleatacsig_rat_to_human_df[i,"Human_geneEnsembl"] <- MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE[MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE$RAT_ENSEMBL_ID %in% ourratens,"HUMAN_ORTHOLOG_ENSEMBL_ID"][1]
    muscleatacsig_rat_to_human_df[i,"Human_geneSymbol"] <- MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE[MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE$RAT_ENSEMBL_ID %in% ourratens,"HUMAN_ORTHOLOG_SYMBOL"][1]
  }
}


muscle_atac_allsigpeak_ORA <- run_ORA(input = unique(muscleatacsig_rat_to_human_df$Human_geneSymbol)[nchar(unique(muscleatacsig_rat_to_human_df$Human_geneSymbol)) > 2],
                                   background = unique(MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE$HUMAN_ORTHOLOG_SYMBOL),
                                   overlap_cutoff = 0.7)

muscle_atac_cluster1_ORA <- run_ORA(input = unique(muscleatacsig_rat_to_human_df[cluster1sigpeaks,"Human_geneSymbol"])[nchar(unique(muscleatacsig_rat_to_human_df[cluster1sigpeaks,"Human_geneSymbol"])) > 2],
                                                                 background = unique(MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE$HUMAN_ORTHOLOG_SYMBOL),
                                                                 overlap_cutoff = 0.7)
muscle_atac_cluster2_ORA <- run_ORA(input = unique(muscleatacsig_rat_to_human_df[cluster2sigpeaks,"Human_geneSymbol"])[nchar(unique(muscleatacsig_rat_to_human_df[cluster2sigpeaks,"Human_geneSymbol"])) > 2],
                                                               background = unique(MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE$HUMAN_ORTHOLOG_SYMBOL),
                                                               overlap_cutoff = 0.7)
muscle_atac_cluster3_ORA <- run_ORA(input = unique(muscleatacsig_rat_to_human_df[cluster3sigpeaks,"Human_geneSymbol"])[nchar(unique(muscleatacsig_rat_to_human_df[cluster3sigpeaks,"Human_geneSymbol"])) > 2],
                                                               background = unique(MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE$HUMAN_ORTHOLOG_SYMBOL),
                                                               overlap_cutoff = 0.7)
muscle_atac_cluster4_ORA <- run_ORA(input = unique(muscleatacsig_rat_to_human_df[cluster4sigpeaks,"Human_geneSymbol"])[nchar(unique(muscleatacsig_rat_to_human_df[cluster4sigpeaks,"Human_geneSymbol"])) > 2],
                                                               background = unique(MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE$HUMAN_ORTHOLOG_SYMBOL),
                                                               overlap_cutoff = 0.7)
muscle_atac_cluster5_ORA <- run_ORA(input = unique(muscleatacsig_rat_to_human_df[cluster5sigpeaks,"Human_geneSymbol"])[nchar(unique(muscleatacsig_rat_to_human_df[cluster5sigpeaks,"Human_geneSymbol"])) > 2],
                                                               background = unique(MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE$HUMAN_ORTHOLOG_SYMBOL),
                                                               overlap_cutoff = 0.7)
rownames(muscle_atac_allsigpeak_ORA) <- muscle_atac_allsigpeak_ORA$set
rownames(muscle_atac_cluster1_ORA) <- muscle_atac_cluster1_ORA$set
rownames(muscle_atac_cluster2_ORA) <- muscle_atac_cluster2_ORA$set
rownames(muscle_atac_cluster3_ORA) <- muscle_atac_cluster3_ORA$set
rownames(muscle_atac_cluster4_ORA) <- muscle_atac_cluster4_ORA$set
rownames(muscle_atac_cluster5_ORA) <- muscle_atac_cluster5_ORA$set

toppathwayset <- unique(c(muscle_atac_allsigpeak_ORA[c(1:10),"set"],
                muscle_atac_cluster1_ORA[c(1:10),"set"],
                muscle_atac_cluster2_ORA[c(1:10),"set"],
                muscle_atac_cluster3_ORA[c(1:10),"set"],
                muscle_atac_cluster4_ORA[c(1:10),"set"],
                muscle_atac_cluster5_ORA[c(1:10),"set"]))

toppathwayORApvalmat <- cbind(muscle_atac_allsigpeak_ORA[toppathwayset,"p_value"],
                              muscle_atac_cluster1_ORA[toppathwayset,"p_value"],
                              muscle_atac_cluster2_ORA[toppathwayset,"p_value"],
                              muscle_atac_cluster3_ORA[toppathwayset,"p_value"],
                              muscle_atac_cluster4_ORA[toppathwayset,"p_value"],
                              muscle_atac_cluster5_ORA[toppathwayset,"p_value"])
colnames(toppathwayORApvalmat) <- c("All Sig Peaks","Cluster1","Cluster2","Cluster3","Cluster4","Cluster5")
rownames(toppathwayORApvalmat) <- toppathwayset

#Figure 2A_P2
png("Figure2A_P2_Muscle ATAC Sig Peak Sets Pathway Enrich_110725.png",width = 12,height = 10,units = "in",res = 600)
pheatmap(-log10(toppathwayORApvalmat),color = colorpanel(101,"white","firebrick"),angle_col = 315)
dev.off()

pdf("Figure2A_P2_Muscle ATAC Sig Peak Sets Pathway Enrich_110725.pdf",width = 12,height = 10)
pheatmap(-log10(toppathwayORApvalmat),color = colorpanel(101,"white","firebrick"),angle_col = 315)
dev.off()


#####
# We initially ran PLIER on the 10k most variable peaks plus the training sig peaks; however,
# we later decided to focus on the training sig peaks alone. We use as prior knowledge the human pathway
# set from the precovid data plus Homer transcription factor motif locations
####

# For the purposes of PLIER, we will first identify the 10k most variable ATAC peaks, and focus our
# PLIER analysis on those 10k, unioned with all the DARs identified in the differential ATAC analysis

gastroatac_sd <- apply(gastro_atacnorm,1,sd)
mostvarpeaks <- union(rownames(gastro_atacnorm)[order(-gastroatac_sd)][1:10000],gastro_atac_training_sig)

peak_by_TF_mat <- matrix(0L,nrow = length(mostvarpeaks),ncol = length(unique(allpeakmotifs$Motif.Name)))
rownames(peak_by_TF_mat) <- mostvarpeaks
colnames(peak_by_TF_mat) <- unique(allpeakmotifs$Motif.Name)

mostvarpeakmotifs <- allpeakmotifs[allpeakmotifs$PositionID %in% mostvarpeaks,]
mostvarpeakmotifs <- mostvarpeakmotifs[mostvarpeakmotifs$MotifScore > 5,]


for(i in 1:dim(mostvarpeakmotifs)[1]){
  if(i%%10000 == 0){
    print(i)
  }
  peak_by_TF_mat[mostvarpeakmotifs[i,"PositionID"],mostvarpeakmotifs[i,"Motif.Name"]] <- 1
}

# First, we must generate the molecular signature pathway matrix to be used as prior knowledge for PLIER
# This uses the MOLECULAR SIGNATURES variable in the MotrpacHumanPreSuspension package.

ourdatabases <- names(MOLECULAR_SIGNATURES)

# For gene-level enrichments, we want the first nine sets of molecular signatures. We will first do this by creating a matrix for each database, then combining together into
# one larger database

biocartagenes <- Reduce(union,MOLECULAR_SIGNATURES$BIOCARTA)
kegg_medicusgenes <- Reduce(union,MOLECULAR_SIGNATURES$KEGG_MEDICUS)
pidgenes <- Reduce(union,MOLECULAR_SIGNATURES$PID)
reactomegenes <- Reduce(union,MOLECULAR_SIGNATURES$REACTOME)
wpgenes <- Reduce(union,MOLECULAR_SIGNATURES$WP)
gobpgenes <- Reduce(union,MOLECULAR_SIGNATURES$GOBP)
goccgenes <- Reduce(union,MOLECULAR_SIGNATURES$GOCC)
gomfgenes <- Reduce(union,MOLECULAR_SIGNATURES$GOMF)
mitocartagenes <- Reduce(union,MOLECULAR_SIGNATURES$MITOCARTA)
combinedgenelist <- Reduce(union,list(biocartagenes,kegg_medicusgenes,pidgenes,reactomegenes,wpgenes,gobpgenes,goccgenes,gomfgenes,mitocartagenes))

biocartamat <- matrix(0L,nrow = length(combinedgenelist),ncol = length(names(MOLECULAR_SIGNATURES$BIOCARTA)))
rownames(biocartamat) <- combinedgenelist
colnames(biocartamat) <- names(MOLECULAR_SIGNATURES$BIOCARTA)
for(i in 1:length(names(MOLECULAR_SIGNATURES$BIOCARTA))){
  biocartamat[MOLECULAR_SIGNATURES$BIOCARTA[i][[1]],i] <- 1
}

kegg_medicusmat <- matrix(0L,nrow = length(combinedgenelist),ncol = length(names(MOLECULAR_SIGNATURES$KEGG_MEDICUS)))
rownames(kegg_medicusmat) <- combinedgenelist
colnames(kegg_medicusmat) <- names(MOLECULAR_SIGNATURES$KEGG_MEDICUS)
for(i in 1:length(names(MOLECULAR_SIGNATURES$KEGG_MEDICUS))){
  kegg_medicusmat[MOLECULAR_SIGNATURES$KEGG_MEDICUS[i][[1]],i] <- 1
}

pidmat <- matrix(0L,nrow = length(combinedgenelist),ncol = length(names(MOLECULAR_SIGNATURES$PID)))
rownames(pidmat) <- combinedgenelist
colnames(pidmat) <- names(MOLECULAR_SIGNATURES$PID)
for(i in 1:length(names(MOLECULAR_SIGNATURES$PID))){
  pidmat[MOLECULAR_SIGNATURES$PID[i][[1]],i] <- 1
}

reactomemat <- matrix(0L,nrow = length(combinedgenelist),ncol = length(names(MOLECULAR_SIGNATURES$REACTOME)))
rownames(reactomemat) <- combinedgenelist
colnames(reactomemat) <- names(MOLECULAR_SIGNATURES$REACTOME)
for(i in 1:length(names(MOLECULAR_SIGNATURES$REACTOME))){
  reactomemat[MOLECULAR_SIGNATURES$REACTOME[i][[1]],i] <- 1
}

wpmat <- matrix(0L,nrow = length(combinedgenelist),ncol = length(names(MOLECULAR_SIGNATURES$WP)))
rownames(wpmat) <- combinedgenelist
colnames(wpmat) <- names(MOLECULAR_SIGNATURES$WP)
for(i in 1:length(names(MOLECULAR_SIGNATURES$WP))){
  wpmat[MOLECULAR_SIGNATURES$WP[i][[1]],i] <- 1
}

gobpmat <- matrix(0L,nrow = length(combinedgenelist),ncol = length(names(MOLECULAR_SIGNATURES$GOBP)))
rownames(gobpmat) <- combinedgenelist
colnames(gobpmat) <- names(MOLECULAR_SIGNATURES$GOBP)
for(i in 1:length(names(MOLECULAR_SIGNATURES$GOBP))){
  gobpmat[MOLECULAR_SIGNATURES$GOBP[i][[1]],i] <- 1
}

goccmat <- matrix(0L,nrow = length(combinedgenelist),ncol = length(names(MOLECULAR_SIGNATURES$GOCC)))
rownames(goccmat) <- combinedgenelist
colnames(goccmat) <- names(MOLECULAR_SIGNATURES$GOCC)
for(i in 1:length(names(MOLECULAR_SIGNATURES$GOCC))){
  goccmat[MOLECULAR_SIGNATURES$GOCC[i][[1]],i] <- 1
}

gomfmat <- matrix(0L,nrow = length(combinedgenelist),ncol = length(names(MOLECULAR_SIGNATURES$GOMF)))
rownames(gomfmat) <- combinedgenelist
colnames(gomfmat) <- names(MOLECULAR_SIGNATURES$GOMF)
for(i in 1:length(names(MOLECULAR_SIGNATURES$GOMF))){
  gomfmat[MOLECULAR_SIGNATURES$GOMF[i][[1]],i] <- 1
}

mitocartamat <- matrix(0L,nrow = length(combinedgenelist),ncol = length(names(MOLECULAR_SIGNATURES$MITOCARTA)))
rownames(mitocartamat) <- combinedgenelist
colnames(mitocartamat) <- names(MOLECULAR_SIGNATURES$MITOCARTA)
for(i in 1:length(names(MOLECULAR_SIGNATURES$MITOCARTA))){
  mitocartamat[MOLECULAR_SIGNATURES$MITOCARTA[i][[1]],i] <- 1
}

finalpathwaymat <- cbind(biocartamat,kegg_medicusmat,pidmat,reactomemat,wpmat,gobpmat,goccmat,gomfmat,mitocartamat)

mol_sig_pathwaymat <- finalpathwaymat


peak_by_pathwaymat <- matrix(0L,nrow = length(mostvarpeaks),ncol = dim(mol_sig_pathwaymat)[2])
rownames(peak_by_pathwaymat) <- mostvarpeaks
colnames(peak_by_pathwaymat) <- colnames(mol_sig_pathwaymat)

for(i in 1:dim(peak_by_pathwaymat)[1]){
  
  if(i%%100 == 0){
    print(i)
  }
  
  ourpeak <- mostvarpeaks[i]
  ourratens <- gastro_atac_all_peakAnnodf[ourpeak,"geneEnsembl"]
  if(ourratens %in% MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE$RAT_ENSEMBL_ID){
    ourhumansym <- MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE[MotrpacRatTraining6moData::RAT_TO_HUMAN_GENE$RAT_ENSEMBL_ID %in% ourratens,"HUMAN_ORTHOLOG_SYMBOL"][1]
    if(ourhumansym %in% rownames(mol_sig_pathwaymat)){
      peak_by_pathwaymat[i,] <- mol_sig_pathwaymat[ourhumansym,]
    }
  }
  
}

mergedpriormat <- cbind(peak_by_TF_mat,peak_by_pathwaymat)

colnames(gastro_atacnorm) <- gsub("X","",colnames(gastro_atacnorm))

gastro_atacmeta <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% colnames(gastro_atacnorm),c("pid","viallabel","sex","group")]
rownames(gastro_atacmeta) <- gastro_atacmeta$viallabel

colnames(gastro_atacmeta) <- c("pid","viallabel","Sex","Group")
gastro_atacmeta$Group <- factor(gastro_atacmeta$Group,levels = c("control","1w","2w","4w","8w"))

gastro_atacmeta <- gastro_atacmeta[order(gastro_atacmeta$Sex,gastro_atacmeta$Group),]
gastro_atacnorm <- gastro_atacnorm[,rownames(gastro_atacmeta)]


pass1b_cols <- list("Group" = MotrpacRatTraining6moData::GROUP_COLORS,
                    "Sex" = MotrpacRatTraining6moData::SEX_COLORS[c("female","male")])

# We trim our feature set farther for PLIER to just be the significant training peaks
mergedpriormat_sig <- mergedpriormat[gastro_atac_training_sig,]
mergedpriormat_sig <- mergedpriormat_sig[,colSums(mergedpriormat_sig) > 0]

gastroatacnorm_sig <- gastro_atacnorm[gastro_atac_training_sig,]

gastroataccombovalue_sig <- num.pc(as.data.frame(gastroatacnorm_sig), seed = 1)
gastroatac_sig.plierResult.all.7=PLIER(as.matrix(gastroatacnorm_sig),mergedpriormat_sig,k=2*gastroataccombovalue_sig,trace=F, frac=0.7, scale=T, seed=1)
saveRDS(gastroatac_sig.plierResult.all.7,file = "gastroatacsig_plierresult110725.rds")

#save.image("pass1b_gastroc_atac_analysis.RData")

for(i in 1:dim(gastroatac_sig.plierResult.all.7$Z)[2]){
  
  top12rows <- rownames(gastroatac_sig.plierResult.all.7$Z[order(-gastroatac_sig.plierResult.all.7$Z[,i]),][1:12,])
  datatoplot <- t(scale(t(gastro_atacnorm[top12rows,])))
  
  png(filename=paste("FigureS2D_gastroatac_sig_PLIER_LV",toString(i),"Plot110725.png",sep=""),width = 8,height = 4,units = "in",res = 600)
  pheatmap(datatoplot,breaks = seq(-5, 5, length.out = 100),color = colorpanel(101,"blue", "white", "red"),show_colnames=F,cluster_cols = FALSE,cluster_rows = FALSE,annotation_col = gastro_atacmeta[,c("Group","Sex")],fontsize=14,annotation_colors = pass1b_cols,cellwidth = 5.5)
  dev.off()
  
  pdf(file=paste("FigureS2D_gastroatac_sig_PLIER_LV",toString(i),"Plot110725.pdf",sep=""),width = 8,height = 4)
  pheatmap(datatoplot,breaks = seq(-5, 5, length.out = 100),color = colorpanel(101,"blue", "white", "red"),show_colnames=F,cluster_cols = FALSE,cluster_rows = FALSE,annotation_col = gastro_atacmeta[,c("Group","Sex")],fontsize=14,annotation_colors = pass1b_cols,cellwidth = 5.5)
  dev.off()
}

for(i in 1:dim(gastroatac_sig.plierResult.all.7$Z)[2]){
  
  ourlv <- paste("LV",i,sep = "")
  metamerge <- gastro_atacmeta
  metamerge$B <- gastroatac_sig.plierResult.all.7$B[i,rownames(gastro_atacmeta)]
  
  ggboxplot(metamerge, x = "Group", y = "B",fill = "Group",facet.by = "Sex") + 
    stat_compare_means(data = metamerge,label = "p.signif",method = "t.test",ref.group = "control",label.y = max(metamerge[,"B"])+(max(metamerge[,"B"]) - min(metamerge[,"B"]))*0.05) + 
    theme(legend.position ="none",axis.text.x=element_text(angle=45, hjust=1)) + 
    xlab("Sample") + ggtitle(paste("LV",toString(i),sep = "")) + ylim(min(metamerge[,"B"])-(max(metamerge[,"B"]) - min(metamerge[,"B"]))*0.1,max(metamerge[,"B"])+(max(metamerge[,"B"]) - min(metamerge[,"B"]))*0.1) + scale_fill_manual(values = pass1b_cols$Group)
  ggsave(filename = paste("Figure2D_muscleatac_sig_boxplot_LV",toString(i),"_082125.png",sep=""),width = 5,height = 3,units = "in",dpi = 600)
  
  ggboxplot(metamerge, x = "Group", y = "B",fill = "Group",facet.by = "Sex") + 
    stat_compare_means(data = metamerge,label = "p.signif",method = "t.test",ref.group = "control",label.y = max(metamerge[,"B"])+(max(metamerge[,"B"]) - min(metamerge[,"B"]))*0.05) + 
    theme(legend.position ="none",axis.text.x=element_text(angle=45, hjust=1)) + 
    xlab("Sample") + ggtitle(paste("LV",toString(i),sep = "")) + ylim(min(metamerge[,"B"])-(max(metamerge[,"B"]) - min(metamerge[,"B"]))*0.1,max(metamerge[,"B"])+(max(metamerge[,"B"]) - min(metamerge[,"B"]))*0.1) + scale_fill_manual(values = pass1b_cols$Group)
  ggsave(filename = paste("Figure2D_muscleatac_sig_boxplot_LV",toString(i),"_082125.pdf",sep=""),width = 5,height = 3)
  
}


#png(file = "Gastro_ATAC_SIG_PLIER_UPLOT.png",width = 8,height = 8,units = "in",res = 600)
#plotU(gastroatac_sig.plierResult.all.7,auc.cutoff = 0.5,angle_col = 315,fdr.cutoff = 0.5)
#dev.off()


# Summary plot of LV response

#gastroatacsigBmat <- gastroatac_sig.plierResult.all.7$B
#rownames(gastroatacsigBmat) <- paste("LV",c(1:dim(gastroatacsigBmat)[1]),sep = "")

#pheatmap(gastroatacsigBmat,cluster_rows = F,cluster_cols = F,breaks = seq(-2, 2, length.out = 100),color = colorpanel(101,"blue", "white", "red"),show_colnames=F,annotation_col = gastro_atacmeta[,c("Group","Sex")],fontsize=14,annotation_colors = pass1b_cols)


# Let's investigate LV10's top markers on Chromosome 7 more closely

toplv10markers <- names(gastroatac_sig.plierResult.all.7$Z[order(-gastroatac_sig.plierResult.all.7$Z[,10]),10][1:12])


######
# Let's incorporate the RNAseq data
####

# Convert ExpressionSet object to DGEList object
dge <- DGEList(
  counts = exprs(TRNSCRPT_GN),
  samples = pData(TRNSCRPT_GN),
  group = pData(TRNSCRPT_GN)$exp_group,
  genes = fData(TRNSCRPT_GN)
)

# log2 counts-per-million reads
cpm_mat <- cpm(dge, log = TRUE)

gastro_rnanorm <- cpm_mat

fw1rnada <- TRNSCRPT_GN_DA$trained_vs_SED[TRNSCRPT_GN_DA$trained_vs_SED$contrast %in% "F_1W - F_SED",]
rownames(fw1rnada) <- fw1rnada$ensembl_gene_id
fw2rnada <- TRNSCRPT_GN_DA$trained_vs_SED[TRNSCRPT_GN_DA$trained_vs_SED$contrast %in% "F_2W - F_SED",]
rownames(fw2rnada) <- fw2rnada$ensembl_gene_id
fw4rnada <- TRNSCRPT_GN_DA$trained_vs_SED[TRNSCRPT_GN_DA$trained_vs_SED$contrast %in% "F_4W - F_SED",]
rownames(fw4rnada) <- fw4rnada$ensembl_gene_id
fw8rnada <- TRNSCRPT_GN_DA$trained_vs_SED[TRNSCRPT_GN_DA$trained_vs_SED$contrast %in% "F_8W - F_SED",]
rownames(fw8rnada) <- fw8rnada$ensembl_gene_id
mw1rnada <- TRNSCRPT_GN_DA$trained_vs_SED[TRNSCRPT_GN_DA$trained_vs_SED$contrast %in% "M_1W - M_SED",]
rownames(mw1rnada) <- mw1rnada$ensembl_gene_id
mw2rnada <- TRNSCRPT_GN_DA$trained_vs_SED[TRNSCRPT_GN_DA$trained_vs_SED$contrast %in% "M_2W - M_SED",]
rownames(mw2rnada) <- mw2rnada$ensembl_gene_id
mw4rnada <- TRNSCRPT_GN_DA$trained_vs_SED[TRNSCRPT_GN_DA$trained_vs_SED$contrast %in% "M_4W - M_SED",]
rownames(mw4rnada) <- mw4rnada$ensembl_gene_id
mw8rnada <- TRNSCRPT_GN_DA$trained_vs_SED[TRNSCRPT_GN_DA$trained_vs_SED$contrast %in% "M_8W - M_SED",]
rownames(mw8rnada) <- mw8rnada$ensembl_gene_id

gastro_rna_l2fc <- cbind(fw1rnada[rownames(gastro_rnanorm),"logFC"],
                         fw2rnada[rownames(gastro_rnanorm),"logFC"],
                         fw4rnada[rownames(gastro_rnanorm),"logFC"],
                         fw8rnada[rownames(gastro_rnanorm),"logFC"],
                         mw1rnada[rownames(gastro_rnanorm),"logFC"],
                         mw2rnada[rownames(gastro_rnanorm),"logFC"],
                         mw4rnada[rownames(gastro_rnanorm),"logFC"],
                         mw8rnada[rownames(gastro_rnanorm),"logFC"])
rownames(gastro_rna_l2fc) <- rownames(gastro_rnanorm)
colnames(gastro_rna_l2fc) <- c("F W1","F W2","F W4","F W8",
                               "M W1","M W2","M W4","M W8")

gastro_rna_pval <- cbind(fw1rnada[rownames(gastro_rnanorm),"adj.P.Val"],
                         fw2rnada[rownames(gastro_rnanorm),"adj.P.Val"],
                         fw4rnada[rownames(gastro_rnanorm),"adj.P.Val"],
                         fw8rnada[rownames(gastro_rnanorm),"adj.P.Val"],
                         mw1rnada[rownames(gastro_rnanorm),"adj.P.Val"],
                         mw2rnada[rownames(gastro_rnanorm),"adj.P.Val"],
                         mw4rnada[rownames(gastro_rnanorm),"adj.P.Val"],
                         mw8rnada[rownames(gastro_rnanorm),"adj.P.Val"])
rownames(gastro_rna_pval) <- rownames(gastro_rnanorm)
colnames(gastro_rna_pval) <- c("F W1","F W2","F W4","F W8",
                               "M W1","M W2","M W4","M W8")

gastro_rna_timewisesig <- rownames(gastro_rna_pval)[apply(gastro_rna_pval,1,min) < 0.05]

#save.image("pass1b_gastroc_atac_analysis.RData")

# We need to display the distribution of adj-p-vals for the genes associated with significant peaks

#hist(-log10(apply(gastro_rna_pval[intersect(gastro_atac_sig_peakAnnodf$geneEnsembl,rownames(gastro_rna_pval)),],1,min)))

# We want to look at correlations between rna matrix and atac matrix. We will convert column
# names to pid's.

gastro_atacnorm_pid <- gastro_atacnorm
gastro_rnanorm_pid <- gastro_rnanorm

for(i in 1:dim(gastro_atacnorm_pid)[2]){
  colnames(gastro_atacnorm_pid)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% colnames(gastro_atacnorm)[i],"pid"]
}

for(i in 1:dim(gastro_rnanorm_pid)[2]){
  colnames(gastro_rnanorm_pid)[i] <- MotrpacRatTraining6moData::PHENO[MotrpacRatTraining6moData::PHENO$viallabel %in% colnames(gastro_rnanorm)[i],"pid"]
}

gastro_atacnorm_pid <- gastro_atacnorm_pid[,intersect(colnames(gastro_atacnorm_pid),colnames(gastro_rnanorm_pid))]
gastro_rnanorm_pid <- gastro_rnanorm_pid[,intersect(colnames(gastro_atacnorm_pid),colnames(gastro_rnanorm_pid))]

gastro_sigatac_cor <- matrix(0L,nrow = length(gastro_atac_training_sig),ncol = 1)
rownames(gastro_sigatac_cor) <- gastro_atac_training_sig
colnames(gastro_sigatac_cor) <- "Correlation"

for(i in 1:length(gastro_atac_training_sig)){
  if(gastro_atac_sig_peakAnnodf[gastro_atac_training_sig[i],"geneEnsembl"] %in% rownames(gastro_rnanorm_pid)){
    gastro_sigatac_cor[i,"Correlation"] <- cor(t(gastro_atacnorm_pid[gastro_atac_training_sig[i],]),gastro_rnanorm_pid[gastro_atac_sig_peakAnnodf[gastro_atac_training_sig[i],"geneEnsembl"],])
  }
}


#save.image("pass1b_gastroc_atac_analysis.RData")

#####
# We want to identify peaks correlated with their associated gene's expression. Let's focus on DEGs first
####

gastro_rna_degs <- rownames(gastro_rna_pval)[apply(gastro_rna_pval,1,min) < 0.05]

gastro_deg_peaks <- gastro_atac_all_peakAnnodf[gastro_atac_all_peakAnnodf$geneEnsembl %in% gastro_rna_degs,]

gastro_deg_dar_pairs <- rownames(gastro_atac_sig_peakAnnodf[gastro_atac_sig_peakAnnodf$geneEnsembl %in% gastro_rna_degs,])

gastro_deg_dar_pair_cor <- data.frame(row.names = gastro_deg_dar_pairs,
                                      "Gene" = gastro_atac_sig_peakAnnodf[gastro_deg_dar_pairs,"geneEnsembl"],
                                      "Symbol" = gastro_atac_sig_peakAnnodf[gastro_deg_dar_pairs,"geneSymbol"])
gastro_deg_dar_pair_cor$Cor <- 0
for(i in 1:dim(gastro_deg_dar_pair_cor)[1]){
  gastro_deg_dar_pair_cor[i,"Cor"] <- cor(t(gastro_atacnorm_pid[rownames(gastro_deg_dar_pair_cor)[i],]),gastro_rnanorm_pid[gastro_atac_sig_peakAnnodf[rownames(gastro_deg_dar_pair_cor)[i],"geneEnsembl"],])
}
gastro_deg_dar_pair_cor$Count <- 0

#Figure 2B_P2
png(file = "Figure2BP2_gastro_deg_dar_pair_cor_boxpot_with_labels_092225.png",width = 8,height = 4,units = "in",res = 600)
ggplot(gastro_deg_dar_pair_cor,mapping = aes(x=Cor,y=Count))+ geom_boxplot(notch = T) + geom_point()  + theme_classic() + theme(text = element_text(size = 20)) + geom_text_repel(data=subset(gastro_deg_dar_pair_cor, Cor > 0.325 | Cor < -0.075),mapping = aes(label=Symbol),angle = 0,max.overlaps = 20,nudge_y = 0.65,size = 6) + xlab("Correlation") + ylab("DAR-DEG Pairs") + scale_y_discrete(labels = NULL, breaks = NULL)
dev.off()

pdf(file = "Figure2BP2_gastro_deg_dar_pair_cor_boxpot_with_labels_092225.pdf",width = 8,height = 4)
ggplot(gastro_deg_dar_pair_cor,mapping = aes(x=Cor,y=Count))+ geom_boxplot(notch = T) + geom_point()  + theme_classic() + theme(text = element_text(size = 20)) + geom_text_repel(data=subset(gastro_deg_dar_pair_cor, Cor > 0.325 | Cor < -0.075),mapping = aes(label=Symbol),angle = 0,max.overlaps = 20,nudge_y = 0.65,size = 6) + xlab("Correlation") + ylab("DAR-DEG Pairs") + scale_y_discrete(labels = NULL, breaks = NULL)
dev.off()

# Scatter plots of Gene-Peak relationships

ann_cols <- list("Tissue" = c("SKM-GN" = "#088c03",
                              "HEART" = "#f28b2f",
                              "HIPPOC" = "#bf7534",
                              "KIDNEY"= "#7553a7",
                              "LIVER" = "#da6c75",
                              "LUNG" = "#04bf8a",
                              "BAT" = "#8c5220",
                              "WAT-SC" = "#214da6"),
                 "Sex" = c("Female" = "#ff6eff",
                           "Male" = "#5555ff"),
                 "Group" = c("control" = "grey",
                             "1w" = "#F7FCB9",
                             "2w" = "#ADDD8E",
                             "4w" = "#238443",
                             "8w" = "#002612"),
                 "Region" = c("3'.UTR" = "#E377C2FF",
                              "5'.UTR" = "#D62728FF",
                              "Distal.Intergenic" = "#BCBD22FF",
                              "Downstream" = "#7F7F7FFF",
                              "Exon" = "#9467BDFF",
                              "Intron" = "#8C564BFF",
                              "Promoter.(1-2kb)" = "#2CA02CFF",
                              "Promoter.(<=1kb)" = "#FF7F0EFF",
                              "Upstream" = "#1F77B4FF"))

gastro_pidmetadf <- data.frame(row.names = colnames(gastro_rnanorm_pid),
                               "pid" = colnames(gastro_rnanorm_pid))
gastro_pidmetadf$Group <- ""
gastro_pidmetadf$Sex <- ""
for(i in 1:dim(gastro_pidmetadf)[1]){
  gastro_pidmetadf[i,"Group"] <- as.character(gastro_atacmeta[gastro_atacmeta$pid %in% gastro_pidmetadf[i,"pid"],"Group"])
  gastro_pidmetadf[i,"Sex"] <- as.character(gastro_atacmeta[gastro_atacmeta$pid %in% gastro_pidmetadf[i,"pid"],"Sex"])
}

ourgenesym <- "Rcan1"
ourpeak <- "chr11:31620793-31621447"
ourgene <- gastro_atac_sig_peakAnnodf[gastro_atac_sig_peakAnnodf$geneSymbol %in% ourgenesym,"geneEnsembl"][1]
ourdf <- data.frame("RNA" = as.vector(t(gastro_rnanorm_pid[ourgene,])),
                    "ATAC" = as.vector(t(gastro_atacnorm_pid[ourpeak,])),
                    "Group" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Group"],
                    "Sex" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Sex"])
png(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".pdf",sep = ""),width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgenesym <- "Il6r"
ourpeak <- "chr2:175346444-175348295"
ourgene <- gastro_atac_sig_peakAnnodf[gastro_atac_sig_peakAnnodf$geneSymbol %in% ourgenesym,"geneEnsembl"][1]
ourdf <- data.frame("RNA" = as.vector(t(gastro_rnanorm_pid[ourgene,])),
                    "ATAC" = as.vector(t(gastro_atacnorm_pid[ourpeak,])),
                    "Group" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Group"],
                    "Sex" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Sex"])
png(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".pdf",sep = ""),width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgenesym <- "Ptp4a3"
ourpeak <- "chr7:105650867-105652826"
ourgene <- gastro_atac_sig_peakAnnodf[gastro_atac_sig_peakAnnodf$geneSymbol %in% ourgenesym,"geneEnsembl"][1]
ourdf <- data.frame("RNA" = as.vector(t(gastro_rnanorm_pid[ourgene,])),
                    "ATAC" = as.vector(t(gastro_atacnorm_pid[ourpeak,])),
                    "Group" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Group"],
                    "Sex" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Sex"])
png(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".pdf",sep = ""),width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgenesym <- "Gpcpd1"
ourpeak <- "chr3:119785567-119786063"
ourgene <- gastro_atac_sig_peakAnnodf[gastro_atac_sig_peakAnnodf$geneSymbol %in% ourgenesym,"geneEnsembl"][1]
ourdf <- data.frame("RNA" = as.vector(t(gastro_rnanorm_pid[ourgene,])),
                    "ATAC" = as.vector(t(gastro_atacnorm_pid[ourpeak,])),
                    "Group" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Group"],
                    "Sex" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Sex"])
png(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".pdf",sep = ""),width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()


#save.image("pass1b_gastroc_atac_analysis.RData")


ourgenesym <- "Fbxo32"
ourpeak <- "chr7:89738198-89738919"
ourgene <- gastro_atac_sig_peakAnnodf[gastro_atac_sig_peakAnnodf$geneSymbol %in% ourgenesym,"geneEnsembl"][1]
ourdf <- data.frame("RNA" = as.vector(t(gastro_rnanorm_pid[ourgene,])),
                    "ATAC" = as.vector(t(gastro_atacnorm_pid[ourpeak,])),
                    "Group" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Group"],
                    "Sex" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Sex"])
png(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".pdf",sep = ""),width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()

ourgenesym <- "Klhl38"
ourpeak <- "chr7:89738198-89738919"
ourgene <- gastro_atac_sig_peakAnnodf[gastro_atac_sig_peakAnnodf$geneSymbol %in% ourgenesym,"geneEnsembl"][1]
ourdf <- data.frame("RNA" = as.vector(t(gastro_rnanorm_pid[ourgene,])),
                    "ATAC" = as.vector(t(gastro_atacnorm_pid[ourpeak,])),
                    "Group" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Group"],
                    "Sex" = gastro_pidmetadf[colnames(gastro_rnanorm_pid),"Sex"])
png(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".png",sep = ""),width = 6,height = 6,units = "in",res = 600)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()
pdf(file = paste("Figure2C_",ourgenesym,"_",gsub(":","-",ourpeak),".pdf",sep = ""),width = 6,height = 6)
ggplot(ourdf,aes(x=RNA,y=ATAC,shape=Sex,color=Group)) + geom_point(size = 5) + theme_classic() + ggtitle(label = paste("SKM-GN ",ourgenesym," vs ",ourpeak,"\nPearson Correlation = ",substr(toString(cor(ourdf$ATAC,ourdf$RNA),sep = ""),start = 1,stop = 6))) + theme(text = element_text(size = 15)) + ylab("ATAC") + xlab("RNA") + scale_color_manual(values = ann_cols$Group)
dev.off()


#i = 10
#ourlv <- paste("LV",i,sep = "")
#metamerge <- gastro_atacmeta
#metamerge$B <- gastroatac_sig.plierResult.all.7$B[i,rownames(gastro_atacmeta)]
#
#ggboxplot(metamerge, x = "Group", y = "B",fill = "Group",facet.by = "Sex") + 
#  stat_compare_means(data = metamerge,label = "p.signif",method = "t.test",ref.group = "control",label.y = max(metamerge[,"B"])+(max(metamerge[,"B"]) - min(metamerge[,"B"]))*0.05) + 
#  theme(legend.position ="none",axis.text.x=element_text(angle=45, hjust=1)) + 
#  xlab("Sample") + ggtitle(paste("LV",toString(i),sep = "")) + ylim(min(metamerge[,"B"])-(max(metamerge[,"B"]) - min(metamerge[,"B"]))*0.1,max(metamerge[,"B"])+(max(metamerge[,"B"]) - min(metamerge[,"B"]))*0.1) + scale_fill_manual(values = pass1b_cols$Group)
#ggsave(filename = paste("muscleatac_sig_boxplot_LV",toString(i),"_shrunk_082125.png",sep=""),width = 8,height = 3,units = "in",dpi = 600)


######
# Making a venn diagram
####

venninput <- list("DAR-Matched Genes" = gastro_atac_sig_peakAnnodf$geneEnsembl,
                  "DEGs" = gastro_rna_degs)
png(file = "Figure2BP1_DAR_vs_DEG_VennDiagram.png",width = 8,height = 4,units = "in",res = 600)
ggvenn(venninput, c("DAR-Matched Genes", "DEGs"),text_size = 8)
dev.off()
pdf(file = "Figure2BP1_DAR_vs_DEG_VennDiagram.pdf",width = 8,height = 4)
ggvenn(venninput, c("DAR-Matched Genes", "DEGs"),text_size = 8)
dev.off()

# If you want to save your progress following running the code
#save.image("pass1b_gastroc_atac_analysis_110725.RData")


#####
# Highlight human acute response for KLHL38 and FBXO32
####

#MotrpacHumanPreSuspension::plot_single_feature(feature = "ENSG00000156804.7",repo_local_dir = "~/GitHub/precovid-analyses/",selected_tissues = "muscle",selected_omes = "transcript-rna-seq",output_file = "FBXO32_rna_plot_093025.png")
#MotrpacHumanPreSuspension::plot_single_feature(feature = "ENSG00000175946.9",repo_local_dir = "~/GitHub/precovid-analyses/",selected_tissues = "muscle",selected_omes = "transcript-rna-seq",output_file = "KLHL38_rna_plot_093025.png")
