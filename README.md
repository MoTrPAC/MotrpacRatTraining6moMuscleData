
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

## Overview

Contains data and code necessary to reproduce analyses that appear in
“Temporal Multi-Omic Analysis Uncovers Sex-Biased Molecular Programs
Underlying Muscle Adaptation to Endurance Training”

Differential analysis of transcriptomics, proteomics, phosphoproteomics,
and metabolomics subcutaneous white adipose tissue -omics datasets.
Enrichment analyses (FGSEA or KSEA) for each -ome. Weighted Gene
Co-expression Network Analysis for all -omes excluding proteomics. The
`data-raw/` folder contains code to reproduce all analyses. Most R
script names match the objects they create, though some scripts create
multiple objects, and so are given more general names. Code to reproduce
figures are found in `vigenttes` for all of the gene set enrichment
analysis, and in the `Other Figures` folders. If you would like to
recreate the entire project from scratch, like in the case of a new
genome annotation, refer to the script `regenerate_figures.sh`, because
some files rely on others in the order of generation.

## Installation

``` r
devtools::install_github("MoTrPAC/MotrpacRatTraining6moMuscleData")
```

## Getting Help

For questions, bug reporting, and data requests for this package, please
<a
href="https://github.com/MoTrPAC/MotrpacRatTraining6moMuscleData/issues"
target="_blank">submit a new issue</a> and include as many details as
possible.

If the concern is related to functions provided in the
<a href="https://github.com/MoTrPAC/MotrpacRatTraining6moMuscleData"
target="_blank">MotrpacRatTraining6moWAT</a> package, please submit an
issue <a
href="https://github.com/MoTrPAC/MotrpacRatTraining6moMuscleData/issues"
target="_blank">here</a> instead.

## Acknowledgements

MoTrPAC is supported by the National Institutes of Health (NIH) Common
Fund through cooperative agreements managed by the National Institute of
Diabetes and Digestive and Kidney Diseases (NIDDK), National Institute
of Arthritis and Musculoskeletal Diseases (NIAMS), and National
Institute on Aging (NIA).

Specifically, the MoTrPAC Study is supported by NIH grants U24OD026629
(Bioinformatics Center), U24DK112349, U24DK112342, U24DK112340,
U24DK112341, U24DK112326, U24DK112331, U24DK112348 (Chemical Analysis
Sites), U01AR071133, U01AR071130, U01AR071124, U01AR071128, U01AR071150,
U01AR071160, U01AR071158 (Clinical Centers), U24AR071113 (Consortium
Coordinating Center), U01AG055133, U01AG055137 and U01AG055135
(PASS/Animal Sites).

## Data Use Agreement

Recipients and their Agents agree that in publications using **any**
data from MoTrPAC public-use data sets they will acknowledge MoTrPAC as
the source of data, including the version number of the data sets used,
e.g.:

- Data used in the preparation of this article were obtained from the
  Molecular Transducers of Physical Activity Consortium (MoTrPAC)
  database, which is available for public access at
  [motrpac-data.org](motrpac-data.org).
- Data used in the preparation of this article were obtained from the
  Molecular Transducers of Physical Activity Consortium (MoTrPAC)
  MotrpacRatTraining6moMuscleData R package \[0.1.1.0\].
