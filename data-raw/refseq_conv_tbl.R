library(Biobase)
library(biomaRt)
library(dplyr)
library(MotrpacRatTraining6moLungData)

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

mart_df1 <- getBM(
  attributes = c(
    "refseq_peptide",
    "refseq_peptide_predicted",
    "uniprot_gn_id",
    "description"
  ),
  mart = mart
)

mart_df2 <- getBM(
  attributes = c(
    "description",
    "rgd_id",
    "rgd_symbol"
  ),
  mart = mart
) %>%
  filter(!grepl("^LOC", rgd_symbol) & !all(grepl("^LOC", rgd_symbol)),
    .by = description
  ) %>%
  arrange(description, rgd_symbol) %>%
  .[!duplicated(.$description), ]

refseq_conv_tbl <- inner_join(mart_df1, mart_df2, by = "description") %>%
  mutate(description = sub(" [[].*[]]$", "", description)) %>%
  distinct() %>%
  mutate(refseq_peptide = ifelse(refseq_peptide == "",
    refseq_peptide_predicted,
    refseq_peptide
  )) %>%
  filter(refseq_peptide != "") %>%
  select(-refseq_peptide_predicted) %>%
  arrange(refseq_peptide, uniprot_gn_id, description) %>%
  group_by(refseq_peptide, description, rgd_id) %>%
  summarise(
    uniprot_gn_id = paste(unique(uniprot_gn_id), collapse = ", "),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  mutate(across(
    .cols = everything(),
    .fns = ~ ifelse(.x == "", NA_character_, .x)
  )) %>%
  filter(!is.na(uniprot_gn_id) & !all(is.na(uniprot_gn_id)),
    .by = refseq_peptide
  ) %>%
  slice_max(
    order_by = description,
    by = c(refseq_peptide, rgd_id, uniprot_gn_id)
  ) %>%
  slice_max(
    order_by = stringr::str_count(uniprot_gn_id, ", "),
    by = refseq_peptide
  )

# Save
write.table(refseq_conv_tbl,
  file = file.path("data-raw", "refseq_conv_tbl.txt"),
  quote = FALSE, row.names = FALSE, sep = "\t"
)
