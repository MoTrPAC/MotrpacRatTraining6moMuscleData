library(MotrpacRatTraining6moMuscleData)
library(dplyr)
library(TMSig)

## C5 GO:BP subcollection (v2023.2.Hs built on Ensembl 110 gene symbols).
## Downloaded from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
human_gene_sets <- file.path(
  "data-raw",
  c(
    "c5.go.v2023.2.Hs.symbols.gmt",
    "c8.all.v2023.2.Hs.symbols.gmt"
  )
) %>%
  lapply(gmt_to_list) %>%
  setNames(c("C5_GO", "C8"))

lengths(human_gene_sets)
# C5_GO    C8
# 10461   830

# We only want the RUBENSTEIN_SKELETAL_MUSCLE sets from the C8 collection
human_gene_sets$C8 <- human_gene_sets$C8 %>%
  .[grepl("^RUBENSTEIN", names(.))]

lengths(human_gene_sets)
# C5_GO    C8
# 10461    11

# Combine sets
names(human_gene_sets) <- NULL
human_gene_sets <- unlist(human_gene_sets, recursive = FALSE)


## Convert from human genes to rat orthologs ----
conv_tbl <- RAT_TO_HUMAN %>%
  select(
    rat_gene_symbol = symbol,
    human_gene_symbol = human_symbol
  ) %>%
  distinct()

rat_gene_sets <- data.frame(
  set = rep(
    names(human_gene_sets),
    lengths(human_gene_sets)
  ),
  human_gene_symbol = unlist(human_gene_sets),
  stringsAsFactors = FALSE
) %>%
  left_join(conv_tbl,
    by = "human_gene_symbol",
    relationship = "many-to-many"
  ) %>%
  filter(!is.na(rat_gene_symbol)) %>%
  distinct() %>%
  filter(length(unique(human_gene_symbol)) >= 5L, .by = set) %>%
  split(x = .$rat_gene_symbol, f = .$set) %>%
  filter_sets(min_size = 5L)

length(rat_gene_sets) # 9293


# Ratio of rat set sizes to their human set sizes
rat_set_sizes <- lengths(rat_gene_sets)
human_set_sizes <- lengths(human_gene_sets)[names(rat_gene_sets)]

ortholog_ratio <- rat_set_sizes / human_set_sizes

hist(ortholog_ratio,
  breaks = seq(0, 1.4, 0.05),
  main = NULL, xlab = "Rat Set Size / Human Set Size"
)


# Save
usethis::use_data(human_gene_sets, overwrite = TRUE, version = 3)
usethis::use_data(rat_gene_sets, overwrite = TRUE, version = 3)
