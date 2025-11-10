#it's worth noting that the names of the bam files is not always
#exactly the same per tissue. I had to make some of the meta files separate
#per tissue. (the 'star_aligned' part)

set_up_SGSeq_meta = function(selected_tissue = NULL){
  rna_meta = MotrpacRatTraining6moData::TRNSCRPT_META
  SGSeq_meta = rna_meta %>%
    dplyr::filter(Tissue == selected_tissue) %>%
    dplyr::select("viallabel", "Seq_length", "Lib_frag_size", "reads") %>%
    dplyr::rename(sample_name = viallabel,
                  read_length = Seq_length,
                  frag_length = Lib_frag_size,
                  lib_size = reads)

  SGSeq_meta$paired_end = TRUE
  SGSeq_meta$sample_name = as.character(SGSeq_meta$sample_name)
  SGSeq_meta$file_bam = paste0(local_bams,"/", SGSeq_meta$sample_name, "_star_aligned.bam") #rn7.2
  return(SGSeq_meta)
}

# I'm using vastus for the numbers that are copy pasted here, just as
# an example
parse_filter_counts = function(sgvc,
                               selected_tissue = NULL)
{
  counts_matrix = as.data.frame(counts(sgvc)) %>%
    mutate(mean_counts = rowMeans(.))
  # summary(counts_matrix$mean_counts)
  # Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
  # 0.00     0.00     0.68    58.88    11.80 81777.74

  data_dictionary = as.data.frame(elementMetadata(sgvc)) #so this actually
  #describes the relevant info for the splicing event
  counts_annotated = cbind(counts_matrix, data_dictionary) %>%
    mutate(variantName = gsub("[/,]", "_", variantName))
  #so unique genes - found in geneName, for vastus: 12375
  #unique variants - found in variantName, for vastus: 54603

  filtered_counts_psi_logit = counts_annotated  %>%
    filter(mean_counts > 5) %>% #arbitrary but less than 5 counts seems very low
    #so we want to find, for each participant(column),
    #percentage of a given isoform (variantName, each row) out of the sum of counts
    #for the variants within a geneName value.Then logit transform
    select(all_of(c(colnames(counts_matrix), "geneName", "variantName"))) %>%
    select(-mean_counts) %>%
    #selecting just to trim the number of columns and info that's not useful
    group_by(geneName) %>%
    mutate(across(where(is.numeric), ~ . / sum(.) * 100)) %>% #this is percentage
    #if there are no counts entirely for the whole gene for some columns, we just force it to 0%.
    mutate(across(where(is.numeric), ~ replace(., . == "NaN", 0))) %>%
    mutate(across(where(is.numeric), ~ log(. / (100 - . )))) %>% #then this is logit of %
    #we cap situations where the logit transform would make things +/- Inf
    mutate(across(where(is.numeric), ~ replace(., . == -Inf, -5))) %>%
    mutate(across(where(is.numeric), ~ replace(., . == Inf, 5))) %>%
    # mutate(count = n()) %>%
    filter(n() > 2 | (n() == 2 & row_number() == 1)) %>%
    #We want to remove situations where there's 1 isoform per gene, because
    #if its 100% spliced in all the time, that's not very useful.
    #then, if there's only 2 variants for an isoform, we remove the first variant
    #because the percentage spliced in has to be 1-the other variant
    ungroup()

  #so now there are:
  # length(unique(filtered_counts_psi_logit$geneName))  [1] 4486
  # length(unique(filtered_counts_psi_logit$variantName))  [1] 13823
  # These are now ready for linear modeling.

  curr_meta = merge(MotrpacRatTraining6moData::TRNSCRPT_META,
                    MotrpacRatTraining6moData::PHENO, by = "viallabel") %>%
    dplyr::filter(viallabel %in% colnames(counts_matrix)) %>%
    dplyr::select(viallabel, sex, tissue, group)

  output_struc = list()
  output_struc[["counts_annotated"]] = counts_annotated
  output_struc[["qc_norm"]] = filtered_counts_psi_logit
  output_struc[["sample_metadata"]] = curr_meta

  return(output_struc)
}

load_parsed_counts = function(selected_tissue = NULL){
  counts_meta_file = readRDS(file.path(here(), "Figures", "Splicing", "data", selected_tissue, "counts_meta.RDS"))
  return(counts_meta_file)
}

load_sgfc = function(selected_tissue = NULL){
  sgfc_ucsc = readRDS(file.path(here(), "Figures", "Splicing", "data", selected_tissue, "sgfc_10psi.RDS"))
  return(sgfc_ucsc)
}
load_sgvc = function(selected_tissue = NULL){
  sgvc = readRDS(file.path(here(), "Figures", "Splicing", "data", selected_tissue, "sgvc_10psi.RDS"))
  return(sgvc)
}

# The Gene mapping is located
# `gs://motrpac-data-hub/pass1b-06/analysis/resources/motrpac-mappings-master_feature_to_gene.txt`

fit_gene_mapping = function(df_with_name){
  master_mapping = read.csv(
    file.path(here(), "Figures", "Splicing", "data",
              "motrpac-mappings-master_feature_to_gene.txt"),
    sep = "\t"
  ) %>%
    distinct(ensembl_gene, gene_symbol, .keep_all = TRUE)
  df_with_name = df_with_name %>%
    mutate(ensembl_gene = sub("_.*", "", variantName)) %>%
    left_join(master_mapping, by = "ensembl_gene") %>%
    mutate(gene_symbol = toupper(gene_symbol)) %>%
    distinct(variantName, .keep_all = TRUE)

  ensembl = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  chromosome_map = getBM(
    attributes = c("ensembl_gene_id", "chromosome_name"),
    filters = "ensembl_gene_id",
    values = unique(df_with_name$ensembl_gene),
    mart = ensembl
  )  %>%
    rename(ensembl_gene = ensembl_gene_id)

  df_annotated = df_with_name %>%
    left_join(chromosome_map, by = "ensembl_gene")

  return(df_annotated)
}

plot_enrichment_bubble = function(camera_enrichment,
                                  num_to_plot = 12,
                                  custom_title = ""){
  plot_df = head(camera_enrichment, num_to_plot) %>%
    mutate(GeneSet = gsub("GOCC|WP|REACTOME|GOMF|GOBP|KEGG", "", GeneSet)) %>%
    mutate(GeneSet = substr(GeneSet, 2, 45))
  p = ggplot(plot_df, aes(x=GeneSet,
                          y=-log10(PValue),
                          size=NGenes,
                          fill=Direction)) +
    geom_point(pch=21, stat="identity") +
    scale_size("Term Size", range=c(3, 8)) +
    coord_flip() +
    ylab(bquote("-log"[10]*"(p-value)")) +
    theme(axis.title.y=element_blank(), panel.grid.major.y=element_line(size=0.25, color="grey")) +
    ggthemes::theme_few() +
    labs(caption = custom_title)
  return(p)
}


load_differential_splicing = function(selected_tissue = c("gastroc", "vastus")){
  lrt_df = lapply(selected_tissue, function(tissue){
    single_tissue_lrt = readRDS(file.path(here(), "Figures", "Splicing", "data", tissue, "LRT_Results.RDS"))
  }) %>% bind_rows()
}
