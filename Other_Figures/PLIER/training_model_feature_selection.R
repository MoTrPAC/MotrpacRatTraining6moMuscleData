library(MotrpacRatTraining6mo)
library(here)
library(Biobase)
library(dplyr)
library(edgeR)


atac_covariates = c("Sample_batch", "peak_enrich.frac_reads_in_peaks.macs2.frip")
trnscrpt_covariates = c("pct_globin", "rin", "pct_umi_dup", "median_5_3_bias")

# ATAC GN data is too large for normal input so it's loaded via load_atac_data()
atac_data = load_atac_data()
transcript_data = MotrpacRatTraining6moMuscleData::TRNSCRPT_GN

atac_results = training_da(
  eset = atac_data,
  tissue = "SKM-GN",
  assay = "ATAC",
  assay_code = "epigen-atac-seq",
  covariates = atac_covariates
)

atac_results = atac_results %>%
  filter(grepl("chr", feature_ID)) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

saveRDS(atac_results, file = file.path(here(),
                                       "Other_Figures",
                                       "PLIER",
                                       "atac_GN_all_peaks.RDS")
)

atac_sig = atac_results %>%
  filter(adj_p_value < 0.1)

saveRDS(atac_sig, file = file.path(here(),
                                   "Other_Figures",
                                   "PLIER",
                                   "atac_GN_sig_peaks.RDS")
)
#792

# oh yea fyi for the 'removed samples' column. There were 2 outliers in RNA but
# these are not included in this table because they were trimmed prior to loading.
# It shouldn't affect anything in terms of DA results here but just for documentation sake I figure I'd mention it

trnscrpt_results = training_da(
  eset = transcript_data,
  tissue = "SKM-GN",
  assay = "TRNSCRPT",
  assay_code = "transcript-rna-seq",
  covariates = trnscrpt_covariates
)

saveRDS(trnscrpt_sig, file = file.path(here(),
                                       "Other_Figures",
                                       "PLIER",
                                       "trsncrpt_GN_all_features.RDS")
)


trnscrpt_results = trnscrpt_results %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))


trnscrpt_sig = trnscrpt_results %>%
  filter(adj_p_value < 0.1)
#586


saveRDS(trnscrpt_sig, file = file.path(here(),
                                       "Other_Figures",
                                       "PLIER",
                                       "trsncrpt_GN_sig_features.RDS")
)
