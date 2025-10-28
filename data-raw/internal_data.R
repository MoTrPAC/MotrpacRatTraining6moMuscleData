library(dplyr)

# Phenotypic data for skeletal muscle
SKM_PHENO <- MotrpacRatTraining6moData::PHENO %>%
  # Bug: redox samples were not assigned a tissue
  filter(grepl("SKM", tissue) | is.na(tissue)) %>%
  select(pid, bid, viallabel, sex, timepoint = group) %>%
  mutate(
    sex = factor(sex,
      levels = c("female", "male"),
      labels = c("Female", "Male")
    ),
    timepoint = ifelse(timepoint == "control",
      "SED", toupper(timepoint)
    ),
    timepoint = factor(timepoint),
    timepoint = relevel(timepoint, ref = "SED")
  ) %>%
  arrange(sex, timepoint) %>%
  # Combine sex and timepoint (easier to specify contrasts with limma)
  mutate(
    exp_group = paste0(substr(sex, 1, 1), "_", timepoint),
    exp_group = factor(exp_group, levels = unique(exp_group))
  ) %>%
  `rownames<-`(.[["viallabel"]])

# Save
usethis::use_data(SKM_PHENO, overwrite = TRUE, version = 3, internal = TRUE)
