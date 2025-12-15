# damon leach
# november 24, 2025

# load libraries
library(here)
library(tidyverse)
library(openxlsx)
library(multcomp)
library(multcompView)
library(patchwork)

# specify file location
here::i_am("Code/ceramide_analysis.R")

# load dataset
cer <- read.xlsx(here("Data","PASS1B_6M_Cer.xlsx"))
names(cer) <- make.unique(names(cer),sep = "__")

# cer_sum <- cer %>%
#   dplyr::select(`F8.ttest`, dplyr::starts_with("F_"), dplyr::starts_with("M_")) %>%
#   dplyr::filter(F8.ttest == "sum") %>%
#   dplyr::select(-F8.ttest)
# 
# hexcer_sum <- cer_sum[2,]
# cer_sum <- cer_sum[1,]
# 
# notoxhexcer_sum <- cer %>%
#   dplyr::select(`F8.ttest`, dplyr::starts_with("F_"), dplyr::starts_with("M_")) %>%
#   dplyr::filter(F8.ttest == "sum not oxidized") %>%
#   dplyr::select(-F8.ttest)

ratio241_240 <- cer %>%
  dplyr::select(`F8.ttest`, dplyr::starts_with("F_"), dplyr::starts_with("M_")) %>%
  dplyr::filter(F8.ttest == "24:1/24:0") %>%
  dplyr::select(-F8.ttest)

# Statistical Analysis 4: Ceramide 24:1/24:0 ratio
ratio_df <- ratio241_240 %>%
  data.frame() %>%
  tibble::remove_rownames() %>%
  t() %>%
  data.frame() %>%
  dplyr::mutate_all(as.numeric) %>%
  tibble::rownames_to_column(var = "Treatment") %>%
  dplyr::rename(sum_val = ".") %>%
  dplyr::mutate(Treatment = stringr::str_split_i(Treatment, "__",1),
                Sex = stringr::str_split_i(Treatment,"_",1),
                TimePt = stringr::str_split_i(Treatment,"_",2),
                Week = stringr::str_remove(TimePt,"T"),
                Week = ifelse(Week == "C",0,Week),
                Week = as.numeric(Week))
ratio_df <- ratio_df %>%
  dplyr::mutate(TimePt_v2 = TimePt) %>%
  dplyr::mutate(TimePt_v2 = stringr::str_remove(TimePt_v2,"T")) %>%
  dplyr::mutate(TimePt_v2 = paste0(TimePt_v2,"W")) %>%
  dplyr::mutate(TimePt_v2 = ifelse(TimePt_v2 == "CW","SED",TimePt_v2)) %>%
  dplyr::mutate(TimePt_v2 = factor(TimePt_v2, levels = c("SED","1W","2W","4W","8W")))

### female
# filter to just female
ratio_df_f <- ratio_df %>%
  dplyr::filter(Sex == "F")

# run anvoa
mod_aov = aov(sum_val ~ TimePt_v2, data = ratio_df_f)
# this tells us that at least one time point is statistically
# different than another
summary(mod_aov)
# more helpful to look at pairwise comparisons
tukey_results <- TukeyHSD(mod_aov)
female_results <- tukey_results
tukey_pvals <- tukey_results$TimePt_v2[,"p adj"]
names(tukey_pvals) <- row.names(tukey_results$TimePt_v2)
tukey_cld <- multcompLetters(tukey_pvals, threshold = 0.05)

tukey_cld_df <- tukey_cld$Letters %>% data.frame() %>%
  tibble::rownames_to_column(var = "TimePt_v2") %>%
  dplyr::rename(CLD = ".") %>%
  dplyr::mutate(TimePt_v2 = factor(TimePt_v2,levels = c("SED","1W","2W","4W","8W")))

p_f <- ratio_df_f %>%
  dplyr::left_join(tukey_cld_df, by = "TimePt_v2") %>%
  ggplot(aes(x = TimePt_v2, y = sum_val, fill = TimePt_v2)) + 
  geom_boxplot() +
  geom_text(aes(label = CLD), y = 0.14) + 
  theme_bw() + 
  ylim(c(0.05,0.15)) +
  labs(x = "",
       y = "24:1/24:0 Ratio",
       title = "Female") + 
  guides(fill = "none") + 
  scale_fill_manual(values = c("SED" = "white",
                               "1W" = "#F7FCB9",
                               "2W" = "#ADDD8E",
                               "4W" = "#238443",
                               "8W" = "#002612"))

### male
ratio_df_m <- ratio_df %>%
  dplyr::filter(Sex == "M")

# run anvoa
mod_aov = aov(sum_val ~ TimePt_v2, data = ratio_df_m)
# this tells us that at least one time point is statistically
# different than another
summary(mod_aov)
# more helpful to look at pairwise comparisons
tukey_results <- TukeyHSD(mod_aov)
male_results <- tukey_results
tukey_pvals <- tukey_results$TimePt_v2[,"p adj"]
names(tukey_pvals) <- row.names(tukey_results$TimePt_v2)
tukey_cld <- multcompLetters(tukey_pvals, threshold = 0.05)

tukey_cld_df <- tukey_cld$Letters %>% data.frame() %>%
  tibble::rownames_to_column(var = "TimePt_v2") %>%
  dplyr::rename(CLD = ".") %>%
  dplyr::mutate(TimePt_v2 = factor(TimePt_v2,levels = c("SED","1W","2W","4W","8W")))

p_m <- ratio_df_m %>%
  dplyr::left_join(tukey_cld_df, by = "TimePt_v2") %>%
  ggplot(aes(x = TimePt_v2, y = sum_val, fill = TimePt_v2)) + 
  geom_boxplot() +
  geom_text(aes(label = CLD), y = 0.14) + 
  theme_bw() + 
  ylim(c(0.05,0.15)) +
  labs(x = "",
       y = "24:1/24:0 Ratio",
       title = "Male") + 
  guides(fill = "none") + 
  scale_fill_manual(values = c("SED" = "white",
                               "1W" = "#F7FCB9",
                               "2W" = "#ADDD8E",
                               "4W" = "#238443",
                               "8W" = "#002612"))

p_both <- p_f|p_m

# save figure
ggsave(filename = here("Figures","ratio_ceramides.pdf"), plot = p_both,dpi = 300,
       height = 5, width = 8)

# file for p-value results
female_res <- female_results$TimePt_v2
male_res <- male_results$TimePt_v2
res_list <- list(female = female_res, male = male_res)
openxlsx::write.xlsx(x = res_list, file = here("Results","ratio_ceramides_statistics.xlsx"))


# the following code is for other possible analyses that are not currently
# being used in the paper

# # Statistical Analysis 1: Ceramide Sum
# cer_df <- cer_sum %>%
#   data.frame() %>%
#   tibble::remove_rownames() %>%
#   t() %>%
#   data.frame() %>%
#   dplyr::mutate_all(as.numeric) %>%
#   tibble::rownames_to_column(var = "Treatment") %>%
#   dplyr::rename(sum_val = ".") %>%
#   dplyr::mutate(Treatment = stringr::str_split_i(Treatment, "__",1),
#                 Sex = stringr::str_split_i(Treatment,"_",1),
#                 TimePt = stringr::str_split_i(Treatment,"_",2),
#                 Week = stringr::str_remove(TimePt,"T"),
#                 Week = ifelse(Week == "C",0,Week),
#                 Week = as.numeric(Week))
# 
# ### female
# # filter to just female
# cer_df_f <- cer_df %>%
#   dplyr::filter(Sex == "F")
# 
# # run anvoa
# mod_aov = aov(sum_val ~ TimePt, data = cer_df_f)
# # this tells us that at least one time point is statistically
# # different than another
# summary(mod_aov)
# # more helpful to look at pairwise comparisons
# tukey_results <- TukeyHSD(mod_aov)
# tukey_pvals <- tukey_results$TimePt[,"p adj"]
# names(tukey_pvals) <- row.names(tukey_results$TimePt)
# tukey_cld <- multcompLetters(tukey_pvals, threshold = 0.05)
# 
# tukey_cld_df <- tukey_cld$Letters %>% data.frame() %>%
#   tibble::rownames_to_column(var = "TimePt") %>%
#   dplyr::rename(CLD = ".") %>%
#   dplyr::mutate(TimePt = factor(TimePt,levels = c("C","T1","T2","T4","T8")))
# 
# p_f <- cer_df_f %>%
#   dplyr::left_join(tukey_cld_df, by = "TimePt") %>%
#   ggplot(aes(x = TimePt, y = sum_val, fill = TimePt)) + 
#   geom_boxplot() +
#   geom_text(aes(label = CLD), y = 5.4) + 
#   theme_bw() + 
#   ylim(c(2.5,5.5)) +
#   labs(x = "",
#        y = "Sum of Ceramides",
#        title = "Female") + 
#   guides(fill = "none") + 
#   scale_fill_manual(values = c("C" = "white",
#                                "T1" = "#F7FCB9",
#                                "T2" = "#ADDD8E",
#                                "T4" = "#238443",
#                                "T8" = "#002612"))
# 
# 
# ### male
# cer_df_m <- cer_df %>%
#   dplyr::filter(Sex == "M")
# 
# # run anvoa
# mod_aov = aov(sum_val ~ TimePt, data = cer_df_m)
# # this tells us that at least one time point is statistically
# # different than another
# summary(mod_aov)
# # more helpful to look at pairwise comparisons
# tukey_results <- TukeyHSD(mod_aov)
# tukey_pvals <- tukey_results$TimePt[,"p adj"]
# names(tukey_pvals) <- row.names(tukey_results$TimePt)
# tukey_cld <- multcompLetters(tukey_pvals, threshold = 0.05)
# 
# tukey_cld_df <- tukey_cld$Letters %>% data.frame() %>%
#   tibble::rownames_to_column(var = "TimePt") %>%
#   dplyr::rename(CLD = ".") %>%
#   dplyr::mutate(TimePt = factor(TimePt,levels = c("C","T1","T2","T4","T8")))
# 
# p_m <- cer_df_m %>%
#   dplyr::left_join(tukey_cld_df, by = "TimePt") %>%
#   ggplot(aes(x = TimePt, y = sum_val, fill = TimePt)) + 
#   geom_boxplot() +
#   geom_text(aes(label = CLD), y = 5.4) + 
#   theme_bw() + 
#   ylim(c(2.5,5.5)) +
#   labs(x = "Time Point",
#        y = "Sum of Ceramides",
#        title = "Male") + 
#   guides(fill = "none") + 
#   scale_fill_manual(values = c("C" = "white",
#                                "T1" = "#F7FCB9",
#                                "T2" = "#ADDD8E",
#                                "T4" = "#238443",
#                                "T8" = "#002612"))
# p_both <- p_f|p_m
# 
# ggsave(filename = here("Figures","sum_ceramides.pdf"), plot = p_both,dpi = 300,
#        height = 5, width = 8)
# 
# # Statistical Analysis 2: HexCeramide Sum
# hexcer_df <- hexcer_sum %>%
#   data.frame() %>%
#   tibble::remove_rownames() %>%
#   t() %>%
#   data.frame() %>%
#   dplyr::mutate_all(as.numeric) %>%
#   tibble::rownames_to_column(var = "Treatment") %>%
#   dplyr::rename(sum_val = ".") %>%
#   dplyr::mutate(Treatment = stringr::str_split_i(Treatment, "__",1),
#                 Sex = stringr::str_split_i(Treatment,"_",1),
#                 TimePt = stringr::str_split_i(Treatment,"_",2),
#                 Week = stringr::str_remove(TimePt,"T"),
#                 Week = ifelse(Week == "C",0,Week),
#                 Week = as.numeric(Week))
# ### female
# # filter to just female
# hexcer_df_f <- hexcer_df %>%
#   dplyr::filter(Sex == "F")
# 
# # run anvoa
# mod_aov = aov(sum_val ~ TimePt, data = hexcer_df_f)
# # this tells us that at least one time point is statistically
# # different than another
# summary(mod_aov)
# # more helpful to look at pairwise comparisons
# tukey_results <- TukeyHSD(mod_aov)
# tukey_pvals <- tukey_results$TimePt[,"p adj"]
# names(tukey_pvals) <- row.names(tukey_results$TimePt)
# tukey_cld <- multcompLetters(tukey_pvals, threshold = 0.05)
# 
# tukey_cld_df <- tukey_cld$Letters %>% data.frame() %>%
#   tibble::rownames_to_column(var = "TimePt") %>%
#   dplyr::rename(CLD = ".") %>%
#   dplyr::mutate(TimePt = factor(TimePt,levels = c("C","T1","T2","T4","T8")))
# 
# p_f <- hexcer_df_f %>%
#   dplyr::left_join(tukey_cld_df, by = "TimePt") %>%
#   ggplot(aes(x = TimePt, y = sum_val, fill = TimePt)) + 
#   geom_boxplot() +
#   geom_text(aes(label = CLD), y = 3.5) + 
#   theme_bw() + 
#   ylim(c(0,4.5)) +
#   labs(x = "Time Point",
#        y = "Sum of HexCeramides",
#        title = "Female") + 
#   guides(fill = "none") + 
#   scale_fill_manual(values = c("C" = "white",
#                                "T1" = "#F7FCB9",
#                                "T2" = "#ADDD8E",
#                                "T4" = "#238443",
#                                "T8" = "#002612"))
# 
# 
# ### male
# hexcer_df_m <- hexcer_df %>%
#   dplyr::filter(Sex == "M")
# 
# # run anvoa
# mod_aov = aov(sum_val ~ TimePt, data = hexcer_df_m)
# # this tells us that at least one time point is statistically
# # different than another
# summary(mod_aov)
# # more helpful to look at pairwise comparisons
# tukey_results <- TukeyHSD(mod_aov)
# tukey_pvals <- tukey_results$TimePt[,"p adj"]
# names(tukey_pvals) <- row.names(tukey_results$TimePt)
# tukey_cld <- multcompLetters(tukey_pvals, threshold = 0.05)
# 
# tukey_cld_df <- tukey_cld$Letters %>% data.frame() %>%
#   tibble::rownames_to_column(var = "TimePt") %>%
#   dplyr::rename(CLD = ".") %>%
#   dplyr::mutate(TimePt = factor(TimePt,levels = c("C","T1","T2","T4","T8")))
# 
# p_m <- hexcer_df_m %>%
#   dplyr::left_join(tukey_cld_df, by = "TimePt") %>%
#   ggplot(aes(x = TimePt, y = sum_val, fill = TimePt)) + 
#   geom_boxplot() +
#   geom_text(aes(label = CLD), y = 3.5) + 
#   theme_bw() + 
#   ylim(c(0,4.5)) +
#   labs(x = "Time Point",
#        y = "Sum of HexCeramides",
#        title = "Male") + 
#   guides(fill = "none") + 
#   scale_fill_manual(values = c("C" = "white",
#                                "T1" = "#F7FCB9",
#                                "T2" = "#ADDD8E",
#                                "T4" = "#238443",
#                                "T8" = "#002612"))
# 
# p_both <- p_f|p_m
# 
# ggsave(filename = here("Figures","sum_hexceramides.pdf"), plot = p_both,dpi = 300,
#        height = 5, width = 8)
# 
# # Statistical Analysis 3: HexCeramide Sum (not oxidized)
# notoxhex_df <- notoxhexcer_sum %>%
#   data.frame() %>%
#   tibble::remove_rownames() %>%
#   t() %>%
#   data.frame() %>%
#   dplyr::mutate_all(as.numeric) %>%
#   tibble::rownames_to_column(var = "Treatment") %>%
#   dplyr::rename(sum_val = ".") %>%
#   dplyr::mutate(Treatment = stringr::str_split_i(Treatment, "__",1),
#                 Sex = stringr::str_split_i(Treatment,"_",1),
#                 TimePt = stringr::str_split_i(Treatment,"_",2),
#                 Week = stringr::str_remove(TimePt,"T"),
#                 Week = ifelse(Week == "C",0,Week),
#                 Week = as.numeric(Week))
# 
# ### female
# # filter to just female
# notoxhex_df_f <- notoxhex_df %>%
#   dplyr::filter(Sex == "F")
# 
# # run anvoa
# mod_aov = aov(sum_val ~ TimePt, data = notoxhex_df_f)
# # this tells us that at least one time point is statistically
# # different than another
# summary(mod_aov)
# # more helpful to look at pairwise comparisons
# tukey_results <- TukeyHSD(mod_aov)
# tukey_pvals <- tukey_results$TimePt[,"p adj"]
# names(tukey_pvals) <- row.names(tukey_results$TimePt)
# tukey_cld <- multcompLetters(tukey_pvals, threshold = 0.05)
# 
# tukey_cld_df <- tukey_cld$Letters %>% data.frame() %>%
#   tibble::rownames_to_column(var = "TimePt") %>%
#   dplyr::rename(CLD = ".") %>%
#   dplyr::mutate(TimePt = factor(TimePt,levels = c("C","T1","T2","T4","T8")))
# 
# p_f <- notoxhex_df_f %>%
#   dplyr::left_join(tukey_cld_df, by = "TimePt") %>%
#   ggplot(aes(x = TimePt, y = sum_val, fill = TimePt)) + 
#   geom_boxplot() +
#   geom_text(aes(label = CLD), y = 2) + 
#   theme_bw() + 
#   ylim(c(0,2.5)) +
#   labs(x = "Time Point",
#        y = "Sum of HexCeramides (Not Oxidized)",
#        title = "Female") + 
#   guides(fill = "none") + 
#   scale_fill_manual(values = c("C" = "white",
#                                "T1" = "#F7FCB9",
#                                "T2" = "#ADDD8E",
#                                "T4" = "#238443",
#                                "T8" = "#002612"))
# 
# 
# ### male
# notoxhex_df_m <- notoxhex_df %>%
#   dplyr::filter(Sex == "M")
# 
# # run anvoa
# mod_aov = aov(sum_val ~ TimePt, data = notoxhex_df_m)
# # this tells us that at least one time point is statistically
# # different than another
# summary(mod_aov)
# # more helpful to look at pairwise comparisons
# tukey_results <- TukeyHSD(mod_aov)
# tukey_pvals <- tukey_results$TimePt[,"p adj"]
# names(tukey_pvals) <- row.names(tukey_results$TimePt)
# tukey_cld <- multcompLetters(tukey_pvals, threshold = 0.05)
# 
# tukey_cld_df <- tukey_cld$Letters %>% data.frame() %>%
#   tibble::rownames_to_column(var = "TimePt") %>%
#   dplyr::rename(CLD = ".") %>%
#   dplyr::mutate(TimePt = factor(TimePt,levels = c("C","T1","T2","T4","T8")))
# 
# p_m <- notoxhex_df_m %>%
#   dplyr::left_join(tukey_cld_df, by = "TimePt") %>%
#   ggplot(aes(x = TimePt, y = sum_val, fill = TimePt)) + 
#   geom_boxplot() +
#   geom_text(aes(label = CLD), y = 2) + 
#   theme_bw() + 
#   ylim(c(0,2.5)) +
#   labs(x = "Time Point",
#        y = "Sum of HexCeramides (Not Oxidized)",
#        title = "Male") + 
#   guides(fill = "none") + 
#   scale_fill_manual(values = c("C" = "white",
#                                "T1" = "#F7FCB9",
#                                "T2" = "#ADDD8E",
#                                "T4" = "#238443",
#                                "T8" = "#002612"))
# 
# p_both <- p_f|p_m
# 
# ggsave(filename = here("Figures","sum_hexceramides_notox.pdf"), plot = p_both,dpi = 300,
#        height = 5, width = 8)
