library(dplyr)

# NOTICE: PSP data is not for commercial use.

# Kinase_Substrate_Dataset obtained from PhosphoSitePlus v6.7.1.1
# https://www.phosphosite.org/staticDownloads.action
# Last Modified Fri Nov 17 08:50:20 EST 2023
psp <- file.path(
  "data-raw",
  "Kinase_Substrate_Dataset.gz"
) %>%
  gzfile() %>%
  read.delim(skip = 2)

human_kinase_sets <- psp %>%
  filter(KIN_ORGANISM == "human", KIN_ORGANISM == SUB_ORGANISM) %>%
  mutate(human_site = paste(SUB_ACC_ID, SUB_MOD_RSD, sep = "-")) %>%
  select(human_site, kinase = GENE) %>%
  distinct() %>%
  # Require at least 5 phosphosites per kinase
  filter(n() >= 3L, .by = kinase) %>%
  unstack()

length(human_kinase_sets) # 319 kinases

hist(lengths(human_kinase_sets)) # most are below 100


# Save
usethis::use_data(human_kinase_sets, overwrite = TRUE, version = 3)
