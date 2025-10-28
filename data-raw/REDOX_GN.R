library(Biobase)
library(dplyr)
library(limma) # removeBatchEffect, plotMDS
library(MotrpacRatTraining6moMuscleData) # normalize_data

refseq_conv_tbl <- file.path("data-raw", "refseq_conv_tbl.txt") %>%
  read.delim()

## Global ratio results ----
global_ratio <- file.path(
  "data-raw",
  "raw_omics_data",
  "MOTRPAC_PASS1B-06_PROT-PR_T55_RN7_Batch2-20240326-results_ratio.txt"
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
  # Remove blank channels
  select(-contains("BLANK", ignore.case = FALSE)) %>%
  `rownames<-`(.[["protein_id"]])


# Expression data
global_expr <- global_ratio %>%
  select(
    starts_with("X90"),
    contains("_TT_", ignore.case = FALSE)
  ) %>%
  as.matrix()

# Feature data
global_fdata <- select(
  global_ratio,
  -all_of(colnames(global_expr))
)

colnames(global_expr) <- sub("^X", "", colnames(global_expr))


## Redox ratio results ----
redox_ratio <- file.path(
  "data-raw",
  "raw_omics_data",
  "MOTRPAC_PASS1B-06_PROT-OX_T55_RN7-20240326-results_ratio.txt"
) %>%
  read.delim() %>%
  mutate(
    refseq_peptide = sub("\\.\\d+$", "", protein_id),
    ptm_id = gsub("c", ";", ptm_id),
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
  # Remove blank channels
  select(-contains("BLANK", ignore.case = FALSE)) %>%
  `rownames<-`(.[["ptm_id"]]) %>%
  # Remove protein from ptm_id
  mutate(ptm_id = sub(".*_", "", ptm_id))

# Expression data
redox_expr <- redox_ratio %>%
  select(
    starts_with("X90"),
    contains("_TT_", ignore.case = FALSE)
  ) %>%
  as.matrix()

# Feature data
redox_fdata <- select(
  redox_ratio,
  -all_of(colnames(redox_expr))
)

colnames(redox_expr) <- sub("^X", "", colnames(redox_expr))

# Order of columns must match
redox_expr <- redox_expr[, colnames(global_expr)]

# Vial metadata
metadata <- file.path(
  "data-raw",
  "raw_omics_data",
  "MOTRPAC_PASS1B-06_T55_OX_PN_20240326_vial_metadata.txt"
) %>%
  read.delim() %>%
  select(viallabel = vial_label, starts_with("tmt"))

p_data <- MotrpacRatTraining6moMuscleData:::SKM_PHENO %>%
  right_join(metadata, by = "viallabel") %>%
  `rownames<-`(.[["viallabel"]]) %>%
  .[colnames(redox_expr), ]


## Create MSnSets ----
# Create ExpressionSet object ----
phenoData <- AnnotatedDataFrame(data = p_data)
global_featureData <- AnnotatedDataFrame(data = global_fdata)
redox_featureData <- AnnotatedDataFrame(data = redox_fdata)

# Global proteomics MSnSet required to scale the redox data
PROT_GN_temp <- ExpressionSet(
  assayData = global_expr,
  phenoData = phenoData,
  featureData = global_featureData
)

REDOX_GN <- ExpressionSet(
  assayData = redox_expr,
  phenoData = phenoData,
  featureData = redox_featureData
)


## Calculate Occupancy ----

# TODO #

# Remove total thiol channels
PROT_GN_temp <- PROT_GN_temp[, !grepl("_TT_", sampleNames(PROT_GN_temp))]
REDOX_GN <- REDOX_GN[, !grepl("_TT_", sampleNames(REDOX_GN))]

## Normalize redox samples by scaled global medians ----
# foo1 <- normalize_redox(m_redox = REDOX_GN,
#                         m_global = PROT_GN_temp,
#                         plexID = "tmt_plex")


## Correct by global + median-MAD normalize
foo1 <- foo2 <- REDOX_GN

exprs(foo1) <- sweep(exprs(foo1),
  MARGIN = 2,
  STATS = apply(exprs(PROT_GN_temp),
    2, median,
    na.rm = TRUE
  ),
  FUN = "-"
)

boxplot(exprs(foo1))

exprs(foo1) <- normalize_data(exprs(foo1), mad = TRUE)
boxplot(exprs(foo1))


exprs(foo2) <- normalize_data(exprs(foo2), mad = TRUE)
boxplot(exprs(foo2))


exprs(REDOX_GN) <- normalize_data(exprs(REDOX_GN), mad = TRUE)


# Correct for TMT plex batch effect ----
plotMDS(exprs(REDOX_GN)) # clear separation by plex

exprs(REDOX_GN) <- removeBatchEffect(
  x = exprs(REDOX_GN),
  batch = REDOX_GN$tmt_plex,
  group = REDOX_GN$exp_group
)


label <- REDOX_GN$timepoint
levels(label) <- c("0", "1", "2", "4", "8")
plot_color <- ifelse(REDOX_GN$sex == "Female", "#ff6eff", "#5555ff")

plotMDS(exprs(REDOX_GN),
  label = label,
  col = plot_color
)
# Samples no longer cluster by plex, and there are no apparent outliers. Some
# timepoint and sex separation.


# Save
usethis::use_data(REDOX_GN, overwrite = TRUE, version = 3)
