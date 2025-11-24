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
library(reshape2)
library(ggpubr)

load("./RD0a_1.1-msnsets_processed_for_phenotype_data.RData")

View(pData(m_ls[["m_glob"]]))

pheno_dat <- data.frame(pData(m_ls[["m_glob"]])) %>% select(-c("pid", "viallabel", "tmt11_channel", "tmt_plex", "exp_group"))

pheno_dat2 <-  melt(pheno_dat, idvar = "bid", timevar = "timepoint", direction = "long")

ggplot(data=pheno_dat2, aes(x=timepoint, y=value, fill = sex)) +
  geom_dotplot(binaxis = "y", 
               stackdir='center', 
               position=position_dodge(0.8)) + 
  theme_classic() +
  stat_summary(fun.data = "mean_sdl", #mean + SD
               fun.args = list(mult = 1),
    position = position_dodge(width = .75)) +
  scale_fill_manual(values = c("magenta", "blue")) +
  theme(legend.key.size = unit(2,"line"), axis.text = element_text(color = "black")) +
  facet_wrap(~factor(variable, levels = unique(pheno_dat2$variable)), scales = "free") 


#look at SED and 8W and only glycogen and citrate_synthase
#switch around parameters to plot by sex rather than timepoint
pheno_dat3 <- pheno_dat2 %>% filter(timepoint == "SED" | timepoint == "8W", variable == "glycogen" | variable == "citrate_synthase"| variable == "delta_vo2max_ml_kg_min")
pheno_dat3 <- pheno_dat3 %>% mutate(grp = paste(timepoint, "_", sex, sep = ""))

#glycogen and CS plots
ggplot(data=pheno_dat3[pheno_dat3$variable == "glycogen" | pheno_dat3$variable == "citrate_synthase", ], 
       aes(x=sex, y=value, fill = timepoint)) +
  geom_dotplot(binaxis = "y", 
               stackdir='center', 
               position=position_dodge(0.8)
               ) + 
  theme_classic() +
  geom_pwc( #https://github.com/kassambara/ggpubr/issues/421 to do p value brackets
    aes(group = timepoint),
    method = "t_test",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) +
  stat_summary(fun.data = "mean_sdl", #mean + SD
               fun.args = list(mult = 1),
               position = position_dodge(width = .75),
               color = "black",
               show.legend = FALSE) +
  scale_fill_manual(values = c("gray", "green3"))  +
  facet_wrap(~factor(variable, levels = unique(pheno_dat2$variable)),
             scales = "free",
             strip.position = "left", #move label strip from top to left
             labeller = as_labeller(c(glycogen = paste("Glycogen (ng/\U00B5L/mg)"), citrate_synthase = paste0("CS (U/\U00B5g protein x 10\U00B3)")))) + #"\U00B5" = mu Greek symbol for micro
  ylab(NULL) +
  guides(fill = guide_legend(title = "Timepoint")) +
  theme(legend.key.size = unit(2.5,"line"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(color = c("magenta", "blue")),
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside")


#look at SED and 8W separately for VO2 Max only (supplemental panel)
ggplot(data=pheno_dat3[pheno_dat3$variable == "delta_vo2max_ml_kg_min", ], aes(x=sex, y=value, fill = timepoint)) +
  geom_dotplot(binaxis = "y", 
               stackdir='center', 
               position=position_dodge(0.8)
  ) + 
  theme_classic() +
  geom_pwc( #https://github.com/kassambara/ggpubr/issues/421 to do p value brackets
    aes(group = timepoint),
    method = "t_test",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) +
  stat_summary(fun.data = "mean_sdl", #mean + SD
               fun.args = list(mult = 1),
               position = position_dodge(width = .75),
               color = "black",
               show.legend = FALSE) +
  scale_fill_manual(values = c("gray", "green3"))  +
  labs(y = bquote(Delta~VO[2]~Max~(ml/kg/min))) +
  guides(fill = guide_legend(title = "Timepoint")) +
  theme(legend.key.size = unit(2.5,"line"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(color = c("magenta", "blue")),
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside")

