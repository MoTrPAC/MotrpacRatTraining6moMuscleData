# Run training (F-test) differential analysis on an ExpressionSet.
# eset:       ExpressionSet whose pData has viallabel, sex (Female/Male), and
#             timepoint (SED/1W/2W/4W/8W) columns — e.g. from load_atac_data()
#             or MotrpacRatTraining6moMuscleData::TRNSCRPT_GN.
# assay/assay_code: strings recorded in the output (e.g. "ATAC"/"epigen-atac-seq"
#             or "TRNSCRPT"/"transcript-rna-seq").
training_da = function(eset,
                       tissue,
                       assay,
                       assay_code,
                       covariates = NULL,
                       outliers = NULL,
                       rdata_outfile = NULL,
                       overwrite = FALSE,
                       verbose = TRUE,
                       n_features = Inf){

  meta_df = as.data.frame(Biobase::pData(eset))
  filt_counts = Biobase::exprs(eset)

  # remove outliers
  removed_outliers = character(0)
  if (!is.null(outliers)) {
    removed_outliers = outliers[outliers %in% rownames(meta_df)]
    meta_df = meta_df[!rownames(meta_df) %in% outliers, ]
    filt_counts = filt_counts[, rownames(meta_df)]
  }

  # derive group column (control/1w/2w/4w/8w) from timepoint (SED/1W/2W/4W/8W)
  meta_df$group = tolower(gsub("SED", "control", as.character(meta_df$timepoint)))

  if("sample_batch" %in% tolower(covariates)){
    which = covariates[grepl("sample_batch", covariates, ignore.case=TRUE)]
    meta_df[,which] = factor(meta_df[,which])
  }
  meta_df$group = factor(meta_df$group, levels=c('control','1w','2w','4w','8w'))
  filt_counts = filt_counts[1:min(nrow(filt_counts),n_features),]

  # split by sex
  sex_res = list()
  ebayes_list = list()
  for(SEX in c('male','female')){

    if(verbose) message(sprintf("Performing F-tests for %s %ss...", tissue, SEX))

    curr_meta = meta_df[tolower(meta_df$sex)==SEX,]
    curr_counts = filt_counts[,curr_meta$viallabel]
    curr_outliers = filter_outliers(tissue=tissue, sex=SEX, outliers=removed_outliers)

    full = paste0(c("~ 1", "group", covariates), collapse=" + ")
    reduced = paste0(c("~ 1", covariates), collapse=" + ")

    # normalize and get voom weights
    design = stats::model.matrix(eval(parse(text=full)), data = curr_meta)
    # check if full rank
    if(!is.fullrank(design)){
      if("sample_batch" %in% tolower(covariates)){
        which = covariates[grepl("sample_batch", covariates, ignore.case=TRUE)]
        curr_cov = covariates[!covariates == which]
        full = paste0(c("~ 1", "group", curr_cov), collapse=" + ")
        reduced = paste0(c("~ 1", curr_cov), collapse=" + ")
        design = model.matrix(eval(parse(text=full)), data = curr_meta)
        warning(sprintf("Sample_batch and group or sex are confounded for %s %s. Removing Sample_batch as a covariate.", tissue, SEX))
      }else{
        stop(sprintf("Model matrix with design %s is not full rank.", full))
      }
    }else{
      curr_cov = covariates
    }
    curr_voom = limma::voom(curr_counts, design=design, normalize.method="quantile")

    limma_model1 = limma::lmFit(curr_voom, design)
    eb_Ftest = limma::eBayes(limma_model1)
    ebayes_list[[SEX]] = eb_Ftest
    res = limma::topTable(eb_Ftest, n=Inf, coef=colnames(design)[grepl('group',colnames(design))], sort.by = "none")
    dt = data.table::data.table(feature_ID=rownames(res),
                                assay=assay,
                                assay_code=assay_code,
                                tissue=tissue,
                                tissue_code=MotrpacRatTraining6moData::TISSUE_ABBREV_TO_CODE[[tissue]],
                                removed_samples=ifelse(length(curr_outliers)>0, paste0(curr_outliers, collapse=','), NA_character_),
                                fscore=res$`F`,
                                p_value = res$P.Value,
                                full_model=gsub(' ','',full),
                                reduced_model=gsub(' ','',reduced))
    sex_res[[SEX]] = dt
  }

  # save to file
  if(!is.null(rdata_outfile)){
    if(overwrite | (!overwrite & !file.exists(rdata_outfile))){
      save(ebayes_list, file=rdata_outfile)
      if(verbose) message(sprintf("'ebayes_list' saved in 'rdata_outfile': %s", rdata_outfile))
    }
  }

  # merge
  merged = data.table::data.table(merge(sex_res[['male']], sex_res[['female']],
                                        by=c("feature_ID","assay","assay_code","tissue","tissue_code"),
                                        suffixes=c('_male','_female')))

  # get a single meta p-value per feature using the male- and female-specific p-values
  merged[,p_value := metap::sumlog(c(p_value_male, p_value_female))$p, by=seq_len(nrow(merged))]

  # reorder columns
  merged = merged[,.(
    feature_ID,
    assay,
    assay_code,
    tissue,
    tissue_code,
    removed_samples_male,
    removed_samples_female,
    fscore_male,
    fscore_female,
    p_value_male,
    p_value_female,
    full_model_male,
    full_model_female,
    reduced_model_male,
    reduced_model_female,
    p_value
  )]

  if(verbose) message("Done.")
  return(as.data.frame(merged))
}


# exact implementation from data-raw/ATAC_GN.R
load_atac_data = function(){

  counts_GN <- file.path(
    here(),
    "data-raw",
    "raw_omics_data",
    "motrpac_pass1b-06_t55-gastrocnemius_epigen-atac-seq_counts_v2.0.txt.gz"
  ) %>%
    read.delim() %>%
    select(-starts_with("X8")) %>%
    mutate(featureName = paste(chrom, start, end, sep = "_")) %>%
    `rownames<-`(.[["featureName"]]) %>%
    select(starts_with("X9"))

  colnames(counts_GN) <- sub("^X", "", colnames(counts_GN))


  ## Sample data ---

  # Library preparation batch and fraction of reads in peaks for ATAC-seq
  covariates <- c("Sample_batch", "peak_enrich.frac_reads_in_peaks.macs2.frip")

  meta <- file.path(
    here(),
    "data-raw",
    "raw_omics_data",
    "motrpac_pass1b-06_epigen-atac-seq_qa-qc-metrics_v2.0.csv"
  ) %>%
    read.delim() %>%
    mutate(vialLabel = as.character(vialLabel)) %>%
    select(vialLabel, all_of(covariates))

  p_data <- MotrpacRatTraining6moMuscleData:::SKM_PHENO %>%
    filter(viallabel %in% colnames(counts_GN)) %>%
    left_join(meta, by = c("viallabel" = "vialLabel")) %>%
    mutate(
      Sample_batch = factor(Sample_batch),
      across(
        .cols = all_of(covariates[2]),
        .fns = ~ scale(.x)[, 1]
      )
    ) %>%
    `rownames<-`(.[["viallabel"]]) %>%
    .[colnames(counts_GN), ]


  # Remove low-count transcripts ----
  # Following 10.12688/f1000research.9005.3
  dge_gn <- DGEList(
    counts = counts_GN,
    samples = p_data,
    group = p_data$exp_group,
    genes = NULL
  )


  # Remove low-count transcripts
  keep <- filterByExpr(dge_gn)
  table(keep)
  #  FALSE   TRUE
  # 239875 856947
  dge_gn <- dge_gn[keep, , keep.lib.sizes = FALSE]

  # Calculate library size normalization factors
  dge_gn <- normLibSizes(dge_gn, method = "TMM")

  plot(dge_gn$samples$norm.factors)

  summary(dge_gn$samples$norm.factors)
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0.7312  0.9125  1.0031  1.0106  1.0939  1.3729


  # Create ExpressionSet object ----
  phenoData <- AnnotatedDataFrame(data = dge_gn$samples)

  ATAC_GN <- ExpressionSet(
    assayData = dge_gn$counts,
    phenoData = phenoData
  )

  return(ATAC_GN)

}
