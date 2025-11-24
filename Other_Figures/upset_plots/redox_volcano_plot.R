library(dplyr)
library(Biobase)
library(openxlsx)
library(MSnSet.utils)
library(MSnID)
library(Biostrings)
library(gtools)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(openxlsx)
library(forcats)
library(MotrpacRatTraining6moMuscleData)

tests <- REDOX_GN_NORM_DA[["trained_vs_SED"]]


#plot p.adj cutoffs only
cutoff <- 0.05
tests2 <- tests %>%
  mutate(point_color = case_when(
    adj.P.Val < cutoff & logFC < 0 ~ "Down", #significantly down
    adj.P.Val < cutoff & logFC > 0 ~ "Up", #significantly up
    TRUE ~ "NS" )) %>% group_by(contrast) %>% arrange(adj.P.Val) %>%
  mutate(P_Adj_Rank_Within_Contrast = row_number(),
         gene_labels = ifelse(P_Adj_Rank_Within_Contrast <=5, gene_symbol, NA)) %>%  ungroup()


plot_counts <- tests2 %>%
  group_by(contrast) %>% count(point_color)
up <- filter(plot_counts, point_color == "Up")
down <- filter(plot_counts, point_color == "Down")


#tests2$contrast <- factor(tests2$contrast, levels = unique(tests$contrast))#reorder factors based on initial test results

#change the order so each row is a timepoint with female and male side by side
tests2$contrast <- factor(tests2$contrast, levels = c("F_1W - F_SED","M_1W - M_SED","F_2W - F_SED","M_2W - M_SED",
                                                      "F_4W - F_SED", "M_4W - M_SED", "F_8W - F_SED", "M_8W - M_SED"))#reorder factors based on initial test results

ggplot(data=tests2, aes(x=logFC, y=-log10(adj.P.Val), col=point_color)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("#5555ff", "red3", "lightgrey"), breaks = c("Down", "Up", "NS")) +
  labs(x = expression(paste("log"[2], "(Fold-Change)"))) +
  facet_wrap(~factor(contrast), nrow = 4) + # plot for each contrast, follow initial stats order
  #add in label to count the significant features 'geom_text' will remove outline
  geom_label(data=down, mapping = aes(x = c(-1.5), y = 4.75, label = paste("n =", n)), color ="#5555ff") +
  geom_label(data=up, mapping = aes(x = c(1.5), y = 4.75, label = paste("n =", n)), color = "red3") +
  geom_hline(yintercept=-log10(cutoff), col="black", lty = "longdash") +
  theme(legend.title = element_blank(), text = element_text(size=15)) +
  scale_y_continuous(limits = c(0,5), sec.axis = sec_axis(~ . * 1, breaks = c(-log10(cutoff)), labels = c(paste(cutoff)))) +
  scale_x_continuous(limits = c(-3,3)) +
  geom_text_repel(data = tests2, aes(label = gene_labels), size =3, box.padding = 0.5, point.padding = 0.2,
                  min.segment.length = 0.1, max.overlaps = 30, show.legend = F) +
  guides(color = guide_legend(override.aes = list(size = 5)))

#export as 7x9 PDF; portrait



