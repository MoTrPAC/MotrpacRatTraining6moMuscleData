library(MotrpacRatTraining6moMuscleData)
library(edgeR)
library(dplyr)
library(Biobase)
library(ggplot2)
library(data.table)


# Matrix of expected counts ----
counts_GN <- file.path(
  "data-raw",
  "raw_omics_data",
  "motrpac_pass1b-06_t55-gastrocnemius_epigen-atac-seq_counts_v2.0.txt.gz"
) %>%
  read.delim() %>%
  select(-starts_with("X8")) %>%
  mutate(featureName = paste(chrom, start, end, sep = "_")) %>%
  `rownames<-`(.[["featureName"]]) %>%
  select(starts_with("X9"))

colnames(counts_GN) <- sub("^X", "", colnames(counts_GN))


## Sample data ---

# Library preparation batch and fraction of reads in peaks for ATAC-seq
covariates <- c("Sample_batch", "peak_enrich.frac_reads_in_peaks.macs2.frip")

meta <- file.path(
  "data-raw",
  "raw_omics_data",
  "motrpac_pass1b-06_epigen-atac-seq_qa-qc-metrics_v2.0.csv"
) %>%
  read.delim() %>%
  mutate(vialLabel = as.character(vialLabel)) %>%
  select(vialLabel, all_of(covariates))

p_data <- MotrpacRatTraining6moMuscleData:::SKM_PHENO %>%
  filter(viallabel %in% colnames(counts_GN)) %>%
  left_join(meta, by = c("viallabel" = "vialLabel")) %>%
  mutate(
    Sample_batch = factor(Sample_batch),
    across(
      .cols = all_of(covariates[2]),
      .fns = ~ scale(.x)[, 1]
    )
  ) %>%
  `rownames<-`(.[["viallabel"]]) %>%
  .[colnames(counts_GN), ]


# Remove low-count transcripts ----
# Following 10.12688/f1000research.9005.3
dge_gn <- DGEList(
  counts = counts_GN,
  samples = p_data,
  group = p_data$exp_group,
  genes = NULL
)

# Remove low-count transcripts
keep <- filterByExpr(dge_gn)
table(keep)
#  FALSE   TRUE
# 239875 856947
dge_gn <- dge_gn[keep, , keep.lib.sizes = FALSE]

# Calculate library size normalization factors
dge_gn <- normLibSizes(dge_gn, method = "TMM")

plot(dge_gn$samples$norm.factors)

summary(dge_gn$samples$norm.factors)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7312  0.9125  1.0031  1.0106  1.0939  1.3729


# Create ExpressionSet object ----
phenoData <- AnnotatedDataFrame(data = dge_gn$samples)

ATAC_GN <- ExpressionSet(
  assayData = dge_gn$counts,
  phenoData = phenoData
)

# dim(ATAC_GN) # 856947 features, 50 samples

# counts <- dge_gn$counts
# lib.size <- getNormLibSizes(dge_gn)
# y <- t(log2(t(counts + 0.5) / (lib.size + 1) * 1e+06))
# y <- normalizeBetweenArrays(y, method = "none")
#
# logCPM <- y
#
# save(logCPM, phenoData,
#      file = "~/Downloads/ATAC-GN_normalized_logCPM.RData",
#      version = 3)

# Check for extreme outliers ----

# # Color by sex
color <- ifelse(ATAC_GN$sex == "Female", "#ff6eff", "#5555ff")
label <- ATAC_GN$timepoint
levels(label) <- c("0", "1", "2", "4", "8")

# Slow
plotMDS(dge_gn, top = 1e5, label = label, col = color)
plotMDS(dge_gn, top = 1e5, label = label, col = color, dim.plot = c(1, 3))
# Potential outliers. Samples appear to separate by sex


# Don't Save -> This data object is 45M, it's quite large but not that crazy. Could be convinced either way. Not going to include for faster times.
# usethis::use_data(ATAC_GN, overwrite = TRUE, version = 3)

## Differential analysis -------------------------------------------------------

# Sex-matched trained vs. sedentary comparisons
timepoint <- rep(c("1W", "2W", "4W", "8W"), times = 2)
sex <- rep(c("F", "M"), each = 4)
contrasts <- sprintf(
  "exp_group%s_%s - exp_group%s_SED",
  sex, timepoint, sex
)

# Include sex comparisons
contrasts <- c(
  contrasts,
  "exp_groupM_SED - exp_groupF_SED",
  "exp_groupM_8W - exp_groupF_8W"
)

contrasts


# Split DA results by contrast type (trained vs. SED, male vs. female) and
# adjust p-values.

split_by_type <- function(x) {
  x %>%
    as.data.frame() %>%
    select(-B) %>%
    filter(!is.na(logFC)) %>%
    # Unlike other results, we remove the featureName column to save memory
    tidyr::separate_wider_delim(
      cols = "featureName", delim = "_",
      names = c("chrom", "start", "end"),
      cols_remove = TRUE
    ) %>%
    mutate(
      across(.cols = c(start, end), as.integer),
      # Simplify contrast names
      contrast = gsub("exp_group", "", contrast),
      # Z-score
      z = zscoreT(t, df = df.total, method = "hill"),
      t = round(t, digits = 3),
      across(
        .cols = c(AveExpr, logFC, df.total),
        .fns = ~ round(.x, digits = 2)
      )
    ) %>%
    # Create indicator column
    mutate(
      group = ifelse(grepl("M.*F", contrast), "MvF",
        "trained_vs_SED"
      ),
      # Convert to factor to preserve order when splitting
      group = factor(group,
        levels = c("trained_vs_SED", "MvF")
      )
    ) %>%
    # Split results into a list of 2 tables
    split(f = .[["group"]]) %>%
    # For each table, adjust p-values across related contrasts
    lapply(function(xi) {
      xi %>%
        select(-group) %>%
        droplevels.data.frame() %>%
        mutate(
          contrast = factor(contrast, levels = unique(contrast)),
          # Adjust p-values across sets of related contrasts
          adj.P.Val = p.adjust(P.Value, method = "BH"),
          across(
            .cols = c(P.Value, adj.P.Val),
            .fns = ~ signif(.x, digits = 2)
          )
        )
    })
}


# Create MArrayLM fitted object
model.str <- paste(c("~ 0", "exp_group", covariates), collapse = " + ")

# Takes ~20-30 minutes
ATAC_GN_FIT <- limmaFit(
  object = ATAC_GN,
  model.str = model.str,
  contrasts = contrasts,
  # Robust prior variance
  trend = FALSE, robust = TRUE,
  # Sample-level quality weights
  var.group = "viallabel",
  plot = FALSE
)

# save(ATAC_GN_FIT, file = "sandbox/ATAC_GN_FIT.RData")

# Save F-test results for trained vs. SED comparisons
ATAC_GN_FTest <- limmaDEA(fit = ATAC_GN_FIT, coef = 1:8, test = "F")
#these files are too large to save as data objects for the package
# saveRDS(ATAC_GN_FTest, file = "~/Downloads/ATAC_GN_DA_Ftest.rds")

# Clean up contrast columns
colnames(ATAC_GN_FTest)[grepl("^coef", colnames(ATAC_GN_FTest))] <-
  gsub(
    "^coef\\.|exp_group", "",
    sub(
      "\\.\\.\\.", " - ",
      colnames(ATAC_GN_FTest)[grepl("^coef", colnames(ATAC_GN_FTest))]
    )
  )


# Differential analysis
ATAC_GN_DA <- limmaDEA(fit = ATAC_GN_FIT) %>%
  split_by_type()
#these files are too large to save as data objects for the package
# saveRDS(ATAC_GN_DA, file = "~/Downloads/ATAC_GN_DA.rds")

# Check p-value histograms
lapply(ATAC_GN_DA, function(x) {
  ggplot(x, aes(x = P.Value)) +
    geom_histogram(breaks = seq(0, 1, 0.05)) +
    facet_wrap(~contrast, ncol = 4) +
    theme_bw()
})
# No issues with p-value histograms.


# Filter based on adjusted p-values --------------------------------------------

## 5% FDR
ATAC_GN_DA_05FDR_norm <- lapply(
  ATAC_GN_DA,
  function(x) filter(x, adj.P.Val < 0.05)
)

vapply(ATAC_GN_DA_05FDR_norm, nrow, integer(1L))
# trained_vs_SED    MvF
#             117   4621

lapply(
  ATAC_GN_DA_05FDR_norm,
  function(y) with(y, table(contrast, adj.P.Val < 0.05))
)

## 10% FDR
ATAC_GN_DA_10FDR_norm <- lapply(
  ATAC_GN_DA,
  function(x) filter(x, adj.P.Val < 0.1)
)

vapply(ATAC_GN_DA_10FDR_norm, nrow, integer(1L))
# trained_vs_SED    MvF
#            362   5602

lapply(
  ATAC_GN_DA_10FDR_norm,
  function(y) with(y, table(contrast, adj.P.Val < 0.1))
)


# # Save
# usethis::use_data(ATAC_GN_DA_05FDR_norm, overwrite = TRUE, version = 3)
# usethis::use_data(ATAC_GN_DA_10FDR_norm, overwrite = TRUE, version = 3)
