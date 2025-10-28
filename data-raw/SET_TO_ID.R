library(MotrpacRatTraining6moMuscleData)
library(dplyr)

# Number of digits for each GeneSetID
n_digits <- floor(log10(length(rat_gene_sets))) + 1L

SET_TO_ID <- data.frame(GeneSet = sort(names(rat_gene_sets))) %>%
  mutate(
    GeneSetID = 1:n(),
    GeneSetID = sprintf(paste0("%0", n_digits, "d"), GeneSetID),
    database = sub("(^[^_]+_).*", "\\1", GeneSet),
    # Shorten descriptions for plotting
    GeneSetShort = sub("^[^_]+_", "", GeneSet),
    GeneSetShort = MotrpacRatTraining6moWAT::cutstr(
      GeneSetShort,
      split = "_", n = 40L
    ),
    GeneSetShort = ifelse(nchar(GeneSet) - nchar(database) >
      nchar(GeneSetShort) + 5L + n_digits,
    sprintf("%s...(%s)", GeneSetShort, GeneSetID),
    sub("^[^_]+_", "", GeneSet)
    ),
    NGenesHuman = lengths(human_gene_sets)[GeneSet]
  ) %>%
  relocate(GeneSetID, .before = everything()) %>%
  select(-database)


# Save
usethis::use_data(SET_TO_ID, overwrite = TRUE, version = 3)
