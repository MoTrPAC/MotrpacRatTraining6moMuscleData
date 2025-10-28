library(Biobase)
library(dplyr)
library(limma) # removeBatchEffect, plotMDS
library(MotrpacRatTraining6moMuscleData) # normalize_data

refseq_conv_tbl <- file.path("data-raw", "refseq_conv_tbl.txt") %>%
  read.delim()

# Ratio results with additional feature annotation columns
ratio_res <- file.path(
  "data-raw",
  "raw_omics_data",
  "MOTRPAC_PASS1B-06_PROT-PR_T55_RN7-20231005_213106-results_ratio.txt"
) %>%
  read.delim() %>%
  mutate(refseq_peptide = sub("\\.\\d+$", "", protein_id)) %>%
  left_join(refseq_conv_tbl, by = "refseq_peptide") %>%
  mutate(prop_nonmissing = rowMeans(
    !is.na(across(.cols = starts_with("X90")))
  )) %>%
  # Remove contaminants, missingness filter (proteins must be detected in at
  # least 2 plexes), only keep RefSeq proteins
  filter(
    prop_nonmissing >= 2 / 6,
    !is_contaminant,
    grepl("^[A-Z]P_", protein_id)
  ) %>%
  select(-is_contaminant, -organism_name, -refseq_peptide) %>%
  relocate(starts_with("X90"), .after = everything()) %>%
  `rownames<-`(.[["protein_id"]])

# Feature data ----
f_data <- select(ratio_res, -starts_with("X90"))

# Expression data ----
expr <- ratio_res %>%
  select(starts_with("X90")) %>%
  as.matrix()
colnames(expr) <- sub("^X", "", colnames(expr))

# Sample data ----
metadata <- file.path(
  "data-raw",
  "raw_omics_data",
  "MOTRPAC_PASS1B-06_T55_PR_PN_20230708_vial_metadata.txt"
) %>%
  read.delim() %>%
  select(viallabel = vial_label, starts_with("tmt"))

p_data <- MotrpacRatTraining6moMuscleData:::SKM_PHENO[colnames(expr), ] %>%
  left_join(metadata, by = "viallabel") %>%
  `rownames<-`(.[["viallabel"]]) %>%
  .[colnames(expr), ]

# Create ExpressionSet object ----
phenoData <- AnnotatedDataFrame(data = p_data)
featureData <- AnnotatedDataFrame(data = f_data)

PROT_GN <- ExpressionSet(
  assayData = expr,
  phenoData = phenoData,
  featureData = featureData
)


# Median-MAD normalization ----
exprs(PROT_GN) <- normalize_data(x = exprs(PROT_GN), mad = TRUE)

# TMT Batch Correction ----
exprs(PROT_GN) <- removeBatchEffect(
  x = exprs(PROT_GN),
  batch = PROT_GN$tmt_plex,
  group = PROT_GN$exp_group
)

label <- PROT_GN$timepoint
levels(label) <- c("0", "1", "2", "4", "8")
label_color <- ifelse(PROT_GN$sex == "Female", "#ff6eff", "#5555ff")

plotMDS(exprs(PROT_GN), labels = label, col = label_color)
# Separation by sex along dim 1. Some timepoint separation.

plotMDS(exprs(PROT_GN),
  labels = label, col = label_color,
  dim.plot = c(1, 3)
)


# Remove outliers ----

# There are 3 samples flagged as outliers due to low loading amounts or poor
# labeling. We will remove these.

outlier_samples <- MotrpacRatTraining6moData::OUTLIERS %>%
  filter(grepl("t55", tissue_code), assay == "PROT") %>%
  pull(viallabel)

PROT_GN <- PROT_GN[, !sampleNames(PROT_GN) %in% outlier_samples]


# Save
usethis::use_data(PROT_GN, overwrite = TRUE, version = 3)
