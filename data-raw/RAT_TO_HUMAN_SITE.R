library(Biobase)
library(dplyr)
library(MotrpacRatTraining6moMuscleData) # PHOSPHO_GN

RAT_TO_HUMAN_SITE <- file.path(
  "data-raw",
  "raw_omics_data",
  "motrpac_pass1b-06_proteomics-ph-RatRN7-to-human-20240228.csv"
) %>%
  read.csv() %>%
  mutate(
    # Replace lower-case amino acids with ';'
    across(
      .cols = everything(),
      .fns = ~ sub(";$", "", gsub("[sty]", ";", .x))
    )
  ) %>%
  filter(ptm_id_rat_refseq %in% rownames(PHOSPHO_GN))


# Save
usethis::use_data(RAT_TO_HUMAN_SITE, overwrite = TRUE, version = 3)
