#!/bin/bash
# ======================================================
# Steps to regenerate all data objects from scratch

# The goal of this document is to record how each data object is created in the proper order, for those
# that are interested in repeating the analysis with their own data objects, or in case
# some aspect of the analysis gets updated (e.g. if someone wants to repeat the analysis for
# when a new genome annotation is built)
#
# For any member of the public just interseted in building the vigenettes, all of the data objects
# are already included in the respository so everything should build as planned. However, if you want to replace
# raw data objects, some of the data objects need to be created in a specific way to support the creation of figures in the way intended by the authors.
#
#
# This document starts assuming that any sequencing level quantification, mass spec quantification is already done. We take those datasets and convert them into a data format suitable for loading.
#
#
# The first thing worth noting is that the analysis for the muscle tissue is done for the vastus and the gastroc, but the gastroc has additional post translational modifications done (e.g. acetyl, phospho, redox) and atac.
# done. Both tissues have trasncript, metab. Some of the ATAC scripts take a while to run, due to larger file sizes,
# but the files should be readily available.
# ======================================================

set -euo pipefail

# Define directories
DATA_RAW="data-raw"
VIGNETTES="vignettes"


echo "------------------------------------------------------"
echo "Step 0: Organizing existing data from the MoTrPAC Consortium"
echo "------------------------------------------------------"

# This is entirely for ease of use, changing names and variables to things that are easier to work with internally, making some

Rscript "data-raw/internal_data.R"

echo "------------------------------------------------------"
echo "Step 1: Building normalized expression objects"
echo "------------------------------------------------------"

# Review the documentation in each of the .R files for rationales of how things
# are done, specifically or the methods text for the general overview of data processing.
#
# Generally speaking, samples considered outliers in the MoTrPAC Nature paper were
# also considered outliers here. There are some other samples that seem, as descriptively visualized
# using an MDS plot, to be outliers. These samples are weighted differently in the linear modeling strategies.

Run all scripts in data-raw to generate normalized omic objects.
for script in $(ls "${DATA_RAW}"/* | grep -E "VL|GN" | grep "\.R$"); do
  echo "Running ${script} ..."
  Rscript "${script}"
done

echo "------------------------------------------------------"
echo "Step 1.1: Building global-protein normalized PTM expression objects"
echo "------------------------------------------------------"

# Rscript -e "rmarkdown::render('${VIGNETTES}/PTM_NORMALIZATION.Rmd', output_format = 'html_document', output_file = NULL, run_pandoc = FALSE)"

echo "------------------------------------------------------"
echo "Step 2: Building gene set and mapping objects"
echo "------------------------------------------------------"

# And then we need to build and annotate gene set objects and mapping
# objects to human. Some of these scripts rely on the normalized expression
# values from part A and/or need to be generated in a specific order.

echo "Running RAT_TO_HUMAN.R ..."
Rscript "data-raw/RAT_TO_HUMAN.R"

echo "Running gene_sets.R ..."
Rscript "data-raw/gene_sets.R"

echo "Running human_kinase_sets.R ..."
Rscript "data-raw/human_kinase_sets.R"

echo "Running RAT_TO_HUMAN_SITE.R ..."
Rscript "data-raw/RAT_TO_HUMAN_SITE.R"

echo "Running SET_TO_ID ..."
Rscript "data-raw/SET_TO_ID.R"


echo "------------------------------------------------------"
echo "Step 3: Running differential analysis"
echo "------------------------------------------------------"

Rscript -e "rmarkdown::render('${VIGNETTES}/differential_analysis.Rmd')"


echo "------------------------------------------------------"
echo "Step 4: Running gene set tests"
echo "------------------------------------------------------"

for script in $(ls "${VIGNETTES}"/* | grep -E 'VL|GN' | grep '\.Rmd$'); do
  echo "Running ${script} ..."

  Rscript -e "rmarkdown::render('${script}')"

done


# Notes for running code for other figures"

# Some of the code for the other figures leverage large files that cant be easily
# uploaded into the github repo, and some other scripts that utilize non-R based
# scripts, but the documentation for the code and where/how to generate upstream files
# is loaded in "Other Figures"

echo "------------------------------------------------------"
echo "Pipeline completed successfully."
echo "------------------------------------------------------"
