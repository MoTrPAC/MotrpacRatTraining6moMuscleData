# Helper functions for human vs MotrPAC comparison analysis

# Load EpiTrain (GSE60590) count data and gene annotations from GEO, merge
# with clinical metadata from the local file EpiTrain_summary_clinical_data.xlsx
# (provided by Malene; must be in the working directory), fit a mixed-effects
# voom model (Timepoint × Sex, blocking on subject), and return topTable DE
# results for the female (women_comp) and male (men_comp) contrasts.
#
# Data sources:
#   - GEO accession GSE60590 (Lindholm et al. 2014, PMID 25373103)
#     Counts:      https://www.ncbi.nlm.nih.gov/geo/download/ acc=GSE60590
#                  file=GSE60590_raw_counts_GRCh38.p13_NCBI.tsv.gz
#     Annotations: https://www.ncbi.nlm.nih.gov/geo/download/ type=rnaseq_counts
#                  file=Human.GRCh38.p13.annot.tsv.gz
#   - EpiTrain_summary_clinical_data.xlsx (local, sheet "Period 1"):
#     BMI and trained leg per subject
#
# Returns a named list:
#   $efit_vfit  — eBayes-fitted MArrayLM object (needed for cameraPR z-scores)
#   $annot      — gene annotation data frame with GeneID and Symbol columns
#   $female     — topTable results for women_comp with Symbol column
#   $male       — topTable results for men_comp with Symbol column
load_epitrain_de <- function() {
  urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"

  gds <- GEOquery::getGEO("GSE60590")
  GSE60590_meta <- Biobase::pData(Biobase::phenoData(gds[[1]]))

  path <- paste(urld, "acc=GSE60590", "file=GSE60590_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&")
  tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

  apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
  annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
  rownames(annot) <- annot$GeneID

  EpiTrain_summary_clinical_data <- readxl::read_excel("EpiTrain_summary_clinical_data.xlsx",
                                                        sheet = "Period 1")
  add_meta <- EpiTrain_summary_clinical_data %>%
    dplyr::select(Subject, BMI, `Trained leg`) %>%
    dplyr::mutate(Subject = as.character(Subject)) %>%
    dplyr::mutate(dominant_leg = dplyr::if_else(stringr::str_detect(`Trained leg`, "\\(D\\)$"),
                                                TRUE, FALSE)) %>%
    dplyr::rename(sample_id = Subject)

  sample_metadata <- GSE60590_meta %>%
    tibble::rownames_to_column("subject") %>%
    dplyr::filter(subject %in% colnames(tbl)) %>%
    dplyr::select(subject, characteristics_ch1.1:characteristics_ch1.5) %>%
    dplyr::rename(gender = characteristics_ch1.1,
                  leg = characteristics_ch1.2,
                  exercise_status = characteristics_ch1.3,
                  timepoint = characteristics_ch1.4,
                  sample_id = characteristics_ch1.5) %>%
    dplyr::mutate(gender = stringr::str_remove(gender, "^gender: "),
                  leg = stringr::str_remove(leg, "^biopsy location: "),
                  exercise_status = stringr::str_remove(exercise_status, "^exercise status: "),
                  timepoint = dplyr::if_else(stringr::str_detect(timepoint, "pre"), "pre", "post"),
                  sample_id = stringr::str_remove(sample_id, "^original sample id: ")) %>%
    dplyr::mutate(sample_id = stringr::str_remove(sample_id, "T\\d{1}")) %>%
    dplyr::arrange(sample_id) %>%
    dplyr::mutate(replicate_id = dplyr::consecutive_id(sample_id)) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(exercise_status = factor(exercise_status, levels = c("untrained", "trained"))) %>%
    dplyr::inner_join(add_meta, dplyr::join_by(sample_id))

  counts <- tbl[, colnames(tbl) %in% sample_metadata$subject]
  counts <- counts[, match(sample_metadata$subject, colnames(counts))]

  x <- edgeR::DGEList(counts=counts, samples=sample_metadata)
  x$samples <- x$samples %>%
    data.frame() %>%
    tidyr::unite(col="gender_ex", c(gender, exercise_status), remove=FALSE)

  keep.exprs <- edgeR::filterByExpr(x, group=x$samples$exercise_status)
  x <- x[keep.exprs,, keep.lib.sizes=FALSE]
  x <- edgeR::calcNormFactors(x, method="TMM")

  design <- model.matrix(~0 + gender_ex + dominant_leg + BMI, x$samples)
  fit_vfit <- edgeR::voomLmFit(x, design, block=x$samples$sample_id,
                                sample.weights=TRUE, plot=FALSE)

  cont.mat <- limma::makeContrasts(
    women_comp = gender_exfemale_trained - gender_exfemale_untrained,
    men_comp = gender_exmale_trained - gender_exmale_untrained,
    interaction_effect = (gender_exfemale_trained - gender_exfemale_untrained) -
                         (gender_exmale_trained - gender_exmale_untrained),
    levels = colnames(design)
  )
  fit_vfit.cont <- limma::contrasts.fit(fit_vfit, cont.mat)
  efit_vfit <- limma::eBayes(fit_vfit.cont)

  extract_results <- function(coef) {
    limma::topTable(efit_vfit, coef=coef, number=nrow(x)) %>%
      tibble::rownames_to_column("GeneID") %>%
      dplyr::mutate(GeneID = as.integer(GeneID)) %>%
      dplyr::left_join(annot[, 1:2], dplyr::join_by(GeneID)) %>%
      dplyr::relocate(Symbol, .after = GeneID)
  }

  list(
    efit_vfit = efit_vfit,
    annot = annot,
    female = extract_results("women_comp"),
    male = extract_results("men_comp")
  )
}

# Load MotrPAC rat gastrocnemius transcriptomics DE results from the package
# dataset TRNSCRPT_GN_DA (generated in vignettes/differential_analysis.Rmd),
# filtering the 8-week trained vs. sedentary contrast for each sex and mapping
# rat gene symbols to uppercase (to match human Symbol convention).
#
# Data source:
#   MotrpacRatTraining6moMuscleData::TRNSCRPT_GN_DA[["trained_vs_SED"]]
#   Contrasts used: "F_8W - F_SED" (female), "M_8W - M_SED" (male)
#
# Returns a named list:
#   $female — DE results for F_8W - F_SED with uppercase Symbol column
#   $male   — DE results for M_8W - M_SED with uppercase Symbol column
load_motrpac_de <- function() {
  data("TRNSCRPT_GN_DA", envir=environment())
  list(
    female = TRNSCRPT_GN_DA[["trained_vs_SED"]] %>%
      dplyr::filter(contrast == "F_8W - F_SED") %>%
      dplyr::mutate(Symbol = toupper(gene_symbol)),
    male = TRNSCRPT_GN_DA[["trained_vs_SED"]] %>%
      dplyr::filter(contrast == "M_8W - M_SED") %>%
      dplyr::mutate(Symbol = toupper(gene_symbol))
  )
}

# Scatter plot with density margins for cameraPR Z-score comparisons.
# data: data frame with a GeneSet column plus x_col and y_col
# top_pathways: optional character vector of GeneSet names to highlight in red
make_zscore_scatter <- function(data, x_col, y_col,
                                 top_pathways = NULL,
                                 title = "", subtitle = "",
                                 x_lab = x_col, y_lab = y_col) {
  p <- ggplot(data, aes(.data[[x_col]], .data[[y_col]], label = GeneSet)) +
    geom_point() +
    geom_abline(linetype = "dashed") +
    theme_bw() +
    labs(title = title, subtitle = subtitle, x = x_lab, y = y_lab) +
    theme(axis.title = element_text(face = "bold", size = 15),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "bold"),
          axis.text = element_text(face = "bold", size = 10))
  if (!is.null(top_pathways))
    p <- p + geom_point(data = data %>% dplyr::filter(GeneSet %in% top_pathways), col = "red")
  p
}
