library(MSnSet.utils)
# devtools::install_github("PNNL-Comp-Mass-Spec/MSnSet.utils")
library(MotrpacRatTrainingPhysiologyData)
#devtools::install_github("MoTrPAC/MotrpacRatTrainingPhysiologyData")
library(MotrpacRatTraining6moMuscleData)
#this package
library(tidyverse)


# Import ------------------------------------------------------------------

# data from Nick
load("data/ACETYL_GN_NORM.rda")
load("data/PHOSPHO_GN_NORM.rda")
load("data/REDOX_GN_NORM.rda")
load("data/PROT_GN.rda")

# meta data
my_package <- "MotrpacRatTrainingPhysiologyData"
data(package = my_package)
name_of_all_datasets <- data.frame(data(package = my_package)$results)$Item
data(list = name_of_all_datasets, package = my_package)


# Process -----------------------------------------------------------------


# convert to msnset
m_acet <- as.MSnSet.ExpressionSet(ACETYL2)
m_phos <- as.MSnSet.ExpressionSet(PHOSPHO2)
m_ox <- as.MSnSet.ExpressionSet(REDOX2)
m_glob <- as.MSnSet.ExpressionSet(PROT_GN)

# check pdata
pData(m_acet)
pData(m_phos)
pData(m_ox)
pData(m_glob)
# some experiment pdata. No metadata.

identical(pData(m_acet), pData(m_phos))
identical(pData(m_acet), pData(m_ox))  # rownames/tmt are different
identical(pData(m_acet), pData(m_ox))
identical(pData(m_acet), pData(m_glob))


# meta data to grab (discussed w/ Gina)
# VO2MAX
# NMR
# ANALYTES
# MUSCLES
# FIBER_TYPES



# VO2MAX
p_vo2max <- VO2MAX %>%
   mutate(delta_vo2max_ml_kg_min = post_vo2max_ml_kg_min - pre_vo2max_ml_kg_min) %>%
   select(pid, delta_vo2max_ml_kg_min)

# NMR
p_nmr <- NMR %>%
   mutate(delta_lean_pct = post_lean_pct - pre_lean_pct,
          delta_fat_pct = post_fat_pct - pre_fat_pct) %>%
   select(pid, delta_lean_pct, delta_fat_pct)

# ANALYTES
p_analytes <- ANALYTES %>%
   select(pid, nefa, glycerol, insulin, corticosterone)

# MUSCLES
p_muscles <- MUSCLES %>%
   select(pid, muscle, glycogen, citrate_synthase, mean_CSA) %>%
   filter(muscle %in% c("LG", "MG")) %>%
   summarise(across(where(is.numeric), \(x){mean(x, na.rm=TRUE)}), .by = "pid")

# FIBER_TYPES
p_fiber_types <- FIBER_TYPES %>%
   select(pid, muscle, type, fiber_area) %>%
   filter(muscle %in% c("LG", "MG")) %>%
   group_by(pid, type) %>%
   summarise(fiber_area_mean = mean(fiber_area, na.rm = TRUE)) %>%
   ungroup() %>%
   mutate(fiber_area_total = sum(fiber_area_mean, na.rm = TRUE), .by = "pid") %>%
   mutate(fiber_area_pct = fiber_area_mean/fiber_area_total,
          type = glue::glue("fiber_type_{type}_area_pct")) %>%
   select(pid, type, fiber_area_pct) %>%
   pivot_wider(names_from = type, values_from = fiber_area_pct)


# combine pdata
pdata <- p_vo2max %>%
   full_join(p_nmr, by = "pid") %>%
   full_join(p_analytes, by = "pid") %>%
   full_join(p_muscles, by = "pid") %>%
   full_join(p_fiber_types, by = "pid") %>%
   mutate(across(where(is.numeric), \(x){if_else(is.nan(x), NA, x)})) %>%
   select(-mean_CSA, -matches("^fiber_type"))

glimpse(pdata)


# add pdata to msnsets
m_ls <- list(m_glob, m_acet, m_phos, m_ox) %>%
   map(\(x){
   pData(x) <- pData(x) %>%
      rownames_to_column() %>%
      mutate(pid = as.character(pid)) %>%
      left_join(pdata, by = "pid") %>%
      mutate(bid = as.character(bid)) %>%
      column_to_rownames()
   return(x)
})
names(m_ls) <- c("m_glob", "m_acet", "m_phos", "m_ox")


# why such poor sample coverage for fiber_types????? 1-4
m_ls[["m_ox"]] %>% pData() %>% View("redox")
m_ls[["m_glob"]] %>% pData() %>% View("glob")

save(m_ls, file = "../output/RD0a_1.1-msnsets_processed.RData")























