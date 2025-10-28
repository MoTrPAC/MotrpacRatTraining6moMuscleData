library(Biobase)
library(dplyr)
library(limma) # plotMDS

# Additional metabolite information, such as RefMet metabolite IDs
refmet_df <- MotrpacRatTraining6moData::METAB_LUNG_DA %>%
  # Remove internal standards
  filter(!grepl("InternalStandard|iSTD", metabolite_refmet)) %>%
  select(
    feature_ID, dataset, site, is_targeted,
    metabolite, metabolite_refmet, mz, rt, neutral_mass
  ) %>%
  mutate(across(
    .cols = c(metabolite, metabolite_refmet),
    .fns = ~ ifelse(is.na(.x), feature_ID, .x)
  )) %>%
  distinct()

# Save metabolite IDs to map to RefMet subclasses
write.table(unique(refmet_df$metabolite_refmet),
  file = file.path("data-raw", "RefMet_metabolite_IDs.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Using
# https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmetF_form.php
# and the "RefMet_metabolite_IDs.txt" file, we obtained additional columns.

# refmet_results.txt generated 2024-03-23

refmet_results <- file.path("data-raw", "refmet_results.txt") %>%
  read.delim() %>%
  rename_all(~ tolower(gsub("\\.", "_", .x))) %>%
  mutate(
    standardized_name = ifelse(standardized_name == "-",
      NA, standardized_name
    ),
    sub_class = ifelse(is.na(standardized_name),
      sub("\\(.*", "", input_name),
      sub_class
    )
  ) %>%
  # formula = neutral mass
  select(-c(formula, exact_mass))

refmet_df <- left_join(refmet_df, refmet_results,
  by = c("metabolite_refmet" = "input_name")
) %>%
  relocate(metabolite, metabolite_refmet, standardized_name,
    .after = feature_ID
  ) %>%
  mutate(across(
    .cols = where(is.character),
    .fns = ~ ifelse(.x == "", NA, .x)
  ))

norm_data <- MotrpacRatTraining6moData::METAB_NORM_DATA_FLAT %>%
  filter(
    tissue == "SKM-GN",
    # Remove internal standards
    !grepl("GTInternalStandard|iSTD", feature_ID)
  ) %>%
  select(-c(feature, assay, tissue)) %>%
  # Remove columns with all missing values
  select(where(~ !all(is.na(.x)))) %>%
  # Include RefMet information and stable metabolite IDs
  left_join(refmet_df, by = c("feature_ID", "dataset")) %>%
  mutate(across(
    .cols = c(metabolite, metabolite_refmet),
    .fns = ~ ifelse(is.na(.x), feature_ID, .x)
  )) %>%
  select(-feature_ID) %>% # same as metabolite column
  relocate(starts_with("1"), .after = everything()) %>%
  as.data.frame() %>%
  `rownames<-`(with(., paste(metabolite_refmet, dataset, sep = "_")))

table(table(norm_data$metabolite_refmet) > 1L)
# FALSE  TRUE
#   932   114

# There are 114 metabolites still detected in multiple platforms.

# group_by(norm_data, metabolite_refmet) %>%
#   filter(n() > 1L) %>%
#   ungroup() %>%
#   arrange(metabolite_refmet) %>%
#   View()

# Split into f_data and sample data
f_data <- select(norm_data, -starts_with("1")) # non-sample columns

norm_data <- select(norm_data, starts_with("1")) %>%
  as.matrix()

p_data <- MotrpacRatTraining6moMuscleData:::SKM_PHENO %>%
  select(-viallabel) %>%
  distinct() %>%
  `rownames<-`(.[["pid"]]) %>%
  .[colnames(norm_data), ]

hist(norm_data, breaks = 30) # approximately Normal
boxplot(norm_data)

# Create ExpressionSet object ----
phenoData <- AnnotatedDataFrame(data = p_data)
featureData <- AnnotatedDataFrame(data = f_data)

METAB_GN <- ExpressionSet(
  assayData = norm_data,
  phenoData = phenoData,
  featureData = featureData
)

color <- ifelse(METAB_GN$sex == "Female", "#ff6eff", "#5555ff")
labels <- METAB_GN$timepoint
levels(labels) <- c("0", "1", "2", "4", "8")
plotMDS(exprs(METAB_GN), labels = labels, col = color)
# Clear sex separation. Some timepoint separation.


# Save
usethis::use_data(METAB_GN, overwrite = TRUE, version = 3)
