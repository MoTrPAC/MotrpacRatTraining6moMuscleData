library(MotrpacRatTraining6moMuscleData)
library(biomaRt)
library(limma) # plotMDS
library(edgeR)
library(dplyr)
library(Biobase)

# Matrix of expected counts ----
counts_GN <- file.path(
  "data-raw",
  "raw_omics_data",
  "motrpac_pass1b-06_t55-gastrocnemius_transcript-rna-seq_rsem-genes-count_v2.0.txt"
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

colnames(counts_GN) <- sub("^X", "", colnames(counts_GN))


## Convert from ensembl gene IDs to gene symbols ----
mart <- useEnsembl(
  biomart = "genes",
  dataset = "rnorvegicus_gene_ensembl",
  version = 110#,
  # mirror = "useast"
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
  values = rownames(counts_GN),
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
# ~11% of the genes begin with "LOC"

prop.table(table(is.na(f_data$gene_symbol)))
# ~4% of the genes are missing


## Sample data ----
# Include additional columns
meta <- MotrpacRatTraining6moData::TRNSCRPT_META %>%
  filter(Tissue == "Gastrocnemius Powder") %>%
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
      if (mean(is.na(cov)) >= 0.05) {
        cov <- NA
      } else {
        # Mean imputation
        cov[is.na(cov)] <- mean(cov)
        # Mean-center and scale
        cov <- scale(cov)[, 1]
      }

      return(cov)
    }
  ))

p_data <- MotrpacRatTraining6moMuscleData:::SKM_PHENO %>%
  filter(viallabel %in% colnames(counts_GN)) %>%
  left_join(meta, by = "viallabel") %>%
  `rownames<-`(.[["viallabel"]])

# Remove internal standard columns
counts_GN <- counts_GN[, rownames(p_data)]
counts_GN = counts_GN[f_data$ensembl_gene_id,]

# Remove low-count transcripts ----
# Following 10.12688/f1000research.9005.3
dge_gn <- DGEList(
  counts = counts_GN,
  samples = p_data,
  group = p_data$exp_group,
  genes = f_data
)

# Remove low-count transcripts
keep <- filterByExpr(dge_gn)
table(keep)
dge_gn <- dge_gn[keep, , keep.lib.sizes = FALSE]

# Calculate library size normalization factors
dge_gn <- normLibSizes(dge_gn, method = "TMM")

plot(dge_gn$samples$norm.factors) # library sizes are approximately equal

# Check for extreme outliers ----

# Color by sex
color <- ifelse(dge_gn$samples$sex == "Female", "#ff6eff", "#5555ff")
# Label by number of weeks of exercise
label <- substr(dge_gn$samples$timepoint, 1, 1)
label <- ifelse(label == "S", "0", label)
plotMDS(dge_gn, top = 1000L, label = label, col = color)
plotMDS(dge_gn, top = 1000L, label = label, col = color, dim.plot = c(1, 3))
# No apparent outliers


# Create ExpressionSet object ----
phenoData <- AnnotatedDataFrame(data = dge_gn$samples)
featureData <- AnnotatedDataFrame(data = dge_gn$genes)

TRNSCRPT_GN <- ExpressionSet(
  assayData = dge_gn$counts,
  phenoData = phenoData,
  featureData = featureData
)

dim(TRNSCRPT_GN) # 13034 features, 50 samples


# Save
usethis::use_data(TRNSCRPT_GN, overwrite = TRUE, version = 3)

#
# counts <- dge_gn$counts
# lib.size <- getNormLibSizes(dge_gn)
# y <- t(log2(t(counts + 0.5) / (lib.size + 1) * 1e+06))
# y <- normalizeBetweenArrays(y, method = "none")
#
# foo <- ExpressionSet(assayData = y,
#                      phenoData = phenoData,
#                      featureData = TRNSCRPT_GN@featureData)
#
# saveRDS(foo,
#      file = "sandbox/TRNSCRPT-GN_normalized_logCPM.rds",
#      version = 3)
