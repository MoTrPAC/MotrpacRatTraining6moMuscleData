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
here::i_am("Code/pc_pe_analysis.R")

# load dataset
pc_pe <- read.xlsx(here("Data","PASS1B_6M_PCPEratio_revised.xlsx"), sheet = "PC_PE")

# analysis
# make the column names unique for data manipulation purposes
names(pc_pe) <- make.unique(names(pc_pe),sep = "__")

# structure the data for analysis
pc_pe_df <- pc_pe[1:2,] %>%
  data.frame() %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "X1") %>%
  t() %>%
  data.frame() %>%
  dplyr::mutate_all(as.numeric) %>%
  tibble::rownames_to_column(var = "Treatment") %>%
  dplyr::mutate(Treatment = stringr::str_split_i(Treatment, "__",1),
                PC_PE_ratio = PC.Sum/PE.sum,
                Sex = stringr::str_split_i(Treatment,"_",1),
                TimePt = stringr::str_split_i(Treatment,"_",2),
                Week = stringr::str_remove(TimePt,"T"),
                Week = ifelse(Week == "C",0,Week),
                Week = as.numeric(Week))

pc_pe_df <- pc_pe_df %>%
  dplyr::mutate(TimePt_v2 = TimePt) %>%
  dplyr::mutate(TimePt_v2 = stringr::str_remove(TimePt_v2,"T")) %>%
  dplyr::mutate(TimePt_v2 = paste0(TimePt_v2,"W")) %>%
  dplyr::mutate(TimePt_v2 = ifelse(TimePt_v2 == "CW","SED",TimePt_v2)) %>%
  dplyr::mutate(TimePt_v2 = factor(TimePt_v2, levels = c("SED","1W","2W","4W","8W")))

# plot of PC/PE by sex
# barplot
# pc_pe_df %>%
#   dplyr::group_by(Treatment, Sex, TimePt) %>%
#   dplyr::summarise(meanRatio = mean(PC_PE_ratio, na.rm = T)) %>%
#   ggplot(aes(x = TimePt, y = meanRatio)) + 
#   geom_col() + 
#   facet_wrap(~Sex) + 
#   theme_bw() + 
#   labs(x = "Time Point",
#        y = "Mean PC/PE Ratio")

# boxplot
# Statistical Analysis
### female
# filter to just female
pc_pe_f <- pc_pe_df %>%
  dplyr::filter(Sex == "F")

# run anvoa
mod_aov = aov(PC_PE_ratio ~ TimePt_v2, data = pc_pe_f)
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

p_f <- pc_pe_f %>%
  dplyr::left_join(tukey_cld_df, by = "TimePt_v2") %>%
  ggplot(aes(x = TimePt_v2, y = PC_PE_ratio, fill = TimePt_v2)) + 
  geom_boxplot() +
  geom_text(aes(label = CLD), y = 6.3) + 
  theme_bw() + 
  ylim(c(3.5,6.5)) +
  labs(x = "",
       y = "PC/PE Ratio",
       title = "Female") + 
  guides(fill = "none") + 
  scale_fill_manual(values = c("SED" = "white",
                               "1W" = "#F7FCB9",
                               "2W" = "#ADDD8E",
                               "4W" = "#238443",
                               "8W" = "#002612"))


### male
pc_pe_m <- pc_pe_df %>%
  dplyr::filter(Sex == "M")

# run anvoa
mod_aov = aov(PC_PE_ratio ~ TimePt_v2, data = pc_pe_m)
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

p_m <- pc_pe_m %>%
  dplyr::left_join(tukey_cld_df, by = "TimePt_v2") %>%
  ggplot(aes(x = TimePt_v2, y = PC_PE_ratio, fill = TimePt_v2)) + 
  geom_boxplot() +
  geom_text(aes(label = CLD), y = 6.3) + 
  theme_bw() + 
  ylim(c(3.5,6.5)) +
  labs(x = "",
       y = "PC/PE Ratio",
       title = "Male") + 
  guides(fill = "none") + 
  scale_fill_manual(values = c("SED" = "white",
                               "1W" = "#F7FCB9",
                               "2W" = "#ADDD8E",
                               "4W" = "#238443",
                               "8W" = "#002612"))


p_both <- p_f|p_m

# save figure
ggsave(filename = here("Figures","pc_pe_ratio.pdf"), plot = p_both,dpi = 300,
       height = 5, width = 8)

# file for p-value results
female_res <- female_results$TimePt_v2
male_res <- male_results$TimePt_v2
res_list <- list(female = female_res, male = male_res)
openxlsx::write.xlsx(x = res_list, file = here("Results","pc_pe_statistics.xlsx"))
