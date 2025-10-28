library(MotrpacRatTraining6moMuscleData)
library(biomaRt)
library(limma) # plotMDS
library(edgeR)
library(dplyr)
library(Biobase)

# Matrix of expected counts ----
counts_VL <- file.path(
  "data-raw",
  "raw_omics_data",
  "motrpac_pass1b-06_t56-vastus-lateralis_transcript-rna-seq_rsem-genes-count_v2.0.txt"
) %>%
  read.delim() %>%
  `rownames<-`(.[["gene_id"]]) %>%
  dplyr::select(-gene_id) %>%
  # Expected counts may have decimal portions
  mutate(across(
    .cols = everything(),
    .fns = ~ floor(as.numeric(.x))
  )) %>%
  as.matrix()

colnames(counts_VL) <- sub("^X", "", colnames(counts_VL))


## Convert from ensembl gene IDs to gene symbols ----
mart <- useEnsembl(
  biomart = "genes",
  dataset = "rnorvegicus_gene_ensembl",
  version = 110
) # mRatBN7.2

# datasets <- listDatasets(mart)
# archives <- listEnsemblArchives()
#
# attributes <- listAttributes(mart)
# filters <- listFilters(mart)

f_data <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "description",
    "external_gene_name"
  ),
  filters = "ensembl_gene_id",
  values = rownames(counts_VL),
  mart = mart
) %>%
  mutate(
    across(
      .cols = where(is.character),
      .fns = ~ ifelse(.x == "", NA, .x)
    ),
    # Extract RGD ID and shorten description
    rgd_id = sub(".*;Acc:(.*)[]]$", "\\1", description),
    rgd_id = as.integer(rgd_id),
    description = sub(" [[].*[]]$", "", description)
  ) %>%
  rename(gene_symbol = external_gene_name) %>%
  `rownames<-`(.[["ensembl_gene_id"]])

prop.table(table(grepl("^LOC", f_data$gene_symbol)))
# ~15% of the genes begin with "LOC"

prop.table(table(is.na(f_data$gene_symbol)))
# ~14% of the genes are missing


## Sample data ----
# Include additional columns
meta <- MotrpacRatTraining6moData::TRNSCRPT_META %>%
  filter(Tissue == "Vastus Lateralis Powder") %>%
  dplyr::select(viallabel,
    rin = RIN, pct_globin,
    pct_umi_dup, median_5_3_bias
  ) %>%
  # Use code from MotrpacRatTraining6mo::fix_covariates
  mutate(across(
    .cols = c(
      rin, pct_globin,
      pct_umi_dup, median_5_3_bias
    ),
    .fns = function(cov) {
      # Mean imputation
      cov[is.na(cov)] <- mean(cov)
      # Mean-center and scale
      cov <- scale(cov)[, 1]

      return(cov)
    }
  ))

p_data <- MotrpacRatTraining6moMuscleData:::SKM_PHENO %>%
  filter(viallabel %in% colnames(counts_VL)) %>%
  left_join(meta, by = "viallabel") %>%
  `rownames<-`(.[["viallabel"]])

# Correct sample order
counts_VL <- counts_VL[, rownames(p_data)]


# Remove low-count transcripts ----
# Following 10.12688/f1000research.9005.3
dge_vl <- DGEList(
  counts = counts_VL,
  samples = p_data,
  group = p_data$exp_group,
  genes = f_data
)

keep <- filterByExpr(dge_vl)
table(keep) # 14176 kept. 16384 discarded.
dge_vl <- dge_vl[keep, , keep.lib.sizes = FALSE]

# Calculate library size normalization factors
dge_vl <- normLibSizes(dge_vl, method = "TMM")

plot(dge_vl$samples$norm.factors)

dge_vl$samples$viallabel[14]
# 90585015603 has a much higher normalization factor, so the library size for
# this sample was much lower than the rest.

# Check for extreme outliers ----

# Color by sex
color <- ifelse(dge_vl$samples$sex == "Female", "#ff6eff", "#5555ff")
# Label by number of weeks of exercise
label <- substr(dge_vl$samples$timepoint, 1, 1)
label <- ifelse(label == "S", "0", label)
plotMDS(dge_vl, top = 1000L, label = label, col = color)
plotMDS(dge_vl, top = 1000L, label = label, col = color, dim.plot = c(1, 3))
# There is an outlying SED female sample. This sample, along with a SED male
# sample, were marked as outliers and removed in the landscape paper. We will do
# the same.

# Remove outliers
outliers <- MotrpacRatTraining6moData::OUTLIERS %>%
  filter(tissue == "SKM-VL") %>%
  pull(viallabel)

dge_vl <- dge_vl[, !dge_vl$samples$viallabel %in% outliers]

# Check MDS plots again
color <- ifelse(dge_vl$samples$sex == "Female", "#ff6eff", "#5555ff")
label <- substr(dge_vl$samples$timepoint, 1, 1)
label <- ifelse(label == "S", "0", label)
plotMDS(dge_vl, top = 1000L, label = label, col = color)
plotMDS(dge_vl, top = 1000L, label = label, col = color, dim.plot = c(1, 3))
# There is an outlying 2W female sample, but we will not remove it. Instead, we
# will include sample-level quality weights in the differential analysis. When
# we do so, this sample gets downweighted.

# Create ExpressionSet object ----
phenoData <- AnnotatedDataFrame(data = dge_vl$samples)
featureData <- AnnotatedDataFrame(data = dge_vl$genes)

TRNSCRPT_VL <- ExpressionSet(
  assayData = dge_vl$counts,
  phenoData = phenoData,
  featureData = featureData
)

dim(TRNSCRPT_VL) # 14176 features, 48 samples


# Save
usethis::use_data(TRNSCRPT_VL,
  internal = FALSE,
  overwrite = TRUE, version = 3
)
