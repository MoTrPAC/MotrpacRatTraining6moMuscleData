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
  "MOTRPAC_PASS1B-06_T55_AC_PN-RN7-lncRNA-20240228_151732-ip-results_ratio.txt"
) %>%
  read.delim() %>%
  mutate(
    refseq_peptide = sub("\\.\\d+$", "", protein_id),
    ptm_id = gsub("k", ";", ptm_id),
    ptm_id = sub("(.*);$", "\\1", ptm_id)
  ) %>%
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
  `rownames<-`(.[["ptm_id"]]) %>%
  # Remove protein from ptm_id
  mutate(ptm_id = sub(".*_", "", ptm_id))


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
  "MOTRPAC_PASS1B-06_T55_AC_PN_20230616_vial_metadata.txt"
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

ACETYL_GN <- ExpressionSet(
  assayData = expr,
  phenoData = phenoData,
  featureData = featureData
)

# Median-MAD normalization ----
exprs(ACETYL_GN) <- normalize_data(x = exprs(ACETYL_GN), mad = TRUE)

# TMT Batch Correction ----
exprs(ACETYL_GN) <- removeBatchEffect(
  x = exprs(ACETYL_GN),
  batch = ACETYL_GN$tmt_plex,
  group = ACETYL_GN$exp_group
)

label <- ACETYL_GN$timepoint
levels(label) <- c("0", "1", "2", "4", "8")
label_color <- ifelse(ACETYL_GN$sex == "Female", "#ff6eff", "#5555ff")

plotMDS(exprs(ACETYL_GN), labels = label, col = label_color)
# No clear separation by sex or timepoint. Several apparently outliers, though
# none were flagged in the landscape paper. These will be downweighted in the
# differential analysis.

plotMDS(exprs(ACETYL_GN),
  labels = label, col = label_color,
  dim.plot = c(1, 3)
)
# No clear separation by sex or timepoint.

# Remove outliers ----

outlier_samples <- MotrpacRatTraining6moData::OUTLIERS %>%
  filter(grepl("t55", tissue_code), assay == "ACETYL") %>%
  pull(viallabel)
# No outliers


# Save
usethis::use_data(ACETYL_GN, overwrite = TRUE, version = 3)
