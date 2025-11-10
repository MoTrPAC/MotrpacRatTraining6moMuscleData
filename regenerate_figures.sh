#!/bin/bash
# ======================================================
# Steps to regenerate all data objects from scratch

The goal of this document is to record how each data object is created in the proper order, for those
that are interested in repeating the analysis with their own data objects, or in case
some aspect of the analysis gets updated (e.g. if someone wants to repeat the analysis for
when a new genome annotation is built)

For any member of the public just interseted in building the vigenettes, all of the data objects
are already included in the respository so everything should build as planned. However, if you want to replace
raw data objects, some of the data objects need to be created in a specific way to support the creation of figures in the way intended by the authors.


This document starts assuming that any sequencing level quantification, mass spec quantification is already done. We take those datasets and convert them into a data format suitable for loading.


The first thing worth noting is that the analysis for the muscle tissue is done for the vastus and the gastroc, but the gastroc has additional post translational modifications done (e.g. acetyl, phospho, redox) and atac.
done. Both tissues have trasncript, metab. Some of the ATAC scripts take a while to run, due to larger file sizes,
but
# ======================================================

set -euo pipefail

# Define directories
DATA_RAW="data-raw"
VIGNETTES="vignettes"

echo "------------------------------------------------------"
echo "Step 1: Building normalized expression objects"
echo "------------------------------------------------------"

Review the documentation in each of the .R files for rationales of how things
are done, specifically or the methods text for the general overview of data processing.

Generally speaking, samples considered outliers in the MoTrPAC Nature paper were
also considered outliers here. There are some other samples that seem, as descriptively visualized
using an MDS plot, to be outliers.

# Run all scripts in data-raw to generate normalized omic objects.
for script in $(ls "${DATA_RAW}"/* | grep -E "VL|GN" | grep "\.R$"); do
  echo "Running ${script} ..."
  Rscript "${script}"
done


These generate all the OME_TISSUE data objects.

For the OME_TISSUE NORM:
What sets have a "NORM" set?
  ls | grep "NORM"
ACETYL_GN_NORM_CAMERA.rda
ACETYL_GN_NORM_DA.rda
ACETYL_GN_NORM_FGSEA.rda
acetyl_GN_NORM_gene_sets.rda
ACETYL_GN_NORM_REDUNDANT_SETS.rda
ACETYL_GN_NORM.rda
PHOSPHO_GN_NORM_CAMERA.rda
PHOSPHO_GN_NORM_FGSEA.rda
phospho_GN_NORM_kinase_sets.rda
PHOSPHO_GN_NORM.rda
redox_GN_NORM_gene_sets.rda
REDOX_GN_NORM_REDUNDANT_SETS.rda
REDOX_GN_NORM.rda

Ok so these "NORM" actually mean that they were normalized for global protein abundance.

Nick Day did this normalization, and Im waiting on some code from him.

echo "------------------------------------------------------"
echo "Step 1.5: Building gene set and mapping objects"
echo "------------------------------------------------------"

And then we need to build and annotate gene set objects and mapping
objects to human. Some of these scripts rely on the normalized expression
values from part A and/or need to be generated in a specific order.


echo "------------------------------------------------------"
echo "Step 2: Running differential analysis"
echo "------------------------------------------------------"

# Run the primary differential analysis vignette
Rscript "${VIGNETTES}/differential_analysis.Rmd"

echo "------------------------------------------------------"
echo "Step 3: Running analyses dependent on differential results"
echo "------------------------------------------------------"

for script in $(ls "${VIGNETTES}"/* | grep -E "VL|GN" | grep "\.Rmd$"); do
  echo "Running ${script} ..."
  Rscript "${script}"
done

echo "------------------------------------------------------"
echo "Step 4: Running fuzzy c-means clustering and gene set tests"
echo "------------------------------------------------------"

# temp -> need to fix this bc i need to merge these Rscripts
#for script in $(ls "${VIGNETTES}"/* | grep -E "fuzzy_cmeans|CMEANS_GST" | grep "\.Rmd$"); do
#  echo "Running ${script} ..."
#  Rscript "${script}"
#done

echo "------------------------------------------------------"
echo "Step 5: Notes for running code for other figures"
echo "------------------------------------------------------"

Some of the code for the other figures leverage large files that cant be easily
uploaded into the github repo, and some other scripts that utilize non-R based
scripts, but the documentation for the code and where/how to generate upstream files
is loaded in "Other Figures"

echo "------------------------------------------------------"
echo "Pipeline completed successfully."
echo "------------------------------------------------------"
