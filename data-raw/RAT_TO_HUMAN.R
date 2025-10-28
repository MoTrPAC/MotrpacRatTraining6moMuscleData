## Convert from rat gene symbols to human gene symbols

library(babelgene) # v22.9
library(Biobase)
library(MotrpacRatTraining6moMuscleData)

# Rat gene symbols from all -omics datasets
rat_genes <- c(
  PROT_GN,
  TRNSCRPT_GN,
  TRNSCRPT_VL,
  ACETYL_GN,
  PHOSPHO_GN,
  REDOX_GN
) %>%
  lapply(function(x) fData(x)$gene_symbol) %>%
  unlist() %>%
  unique()

rat_genes <- rat_genes[!is.na(rat_genes)]
length(rat_genes) # 14104

# Rat to human gene conversion table
RAT_TO_HUMAN <- orthologs(
  genes = rat_genes,
  species = "Rattus norvegicus",
  human = FALSE
)
dim(RAT_TO_HUMAN) # 13023   9

hist(table(RAT_TO_HUMAN$human_symbol))
prop.table(table(table(RAT_TO_HUMAN$human_symbol) > 1))
#      ALSE        TRUE
# 0.98721646 0.01278354


# Save
usethis::use_data(RAT_TO_HUMAN, overwrite = TRUE, version = 3)
