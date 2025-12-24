## ExpressionSet Objects =======================================================

#' @title Skeletal muscle ExpressionSet objects
#'
#' @description Object of class
#'   \code{\link[Biobase:ExpressionSet-class]{Biobase::ExpressionSet}}
#'   containing -omics data for differential analysis.
#'
#' @details The proteomics expression datasets contain log\eqn{_2}-transformed
#'   TMT intensities relative to a global reference in each plex. These values
#'   are median-MAD-normalized and TMT-plex-batch-corrected with
#'   \code{\link[limma:removeBatchEffect]{limma::removeBatchEffect}}. Each of
#'   the datasets tagged `NORM` reflect PTM adjustment by global proteomics.
#'
#' @keywords datasets
#'
#' @name SKM_EXP

#' @rdname SKM_EXP
#' @format NULL
#' @usage NULL
"TRNSCRPT_GN"

#' @rdname SKM_EXP
#' @format NULL
#' @usage NULL
"PROT_GN"

#' @rdname SKM_EXP
#' @format NULL
#' @usage NULL
"PHOSPHO_GN"

#' @rdname SKM_EXP
#' @format NULL
#' @usage NULL
"PHOSPHO_GN_NORM"

#' @rdname SKM_EXP
#' @format NULL
#' @usage NULL
"ACETYL_GN"

#' @rdname SKM_EXP
#' @format NULL
#' @usage NULL
"ACETYL_GN_NORM"

#' @rdname SKM_EXP
#' @format NULL
#' @usage NULL
"REDOX_GN"

#' @rdname SKM_EXP
#' @format NULL
#' @usage NULL
"REDOX_GN_NORM"

#' @rdname SKM_EXP
#' @format NULL
#' @usage NULL
"METAB_GN"

#' @rdname SKM_EXP
#' @format NULL
#' @usage NULL
"TRNSCRPT_VL"

#' @rdname SKM_EXP
#' @format NULL
#' @usage NULL
"METAB_VL"


## Differential Analysis =======================================================

#' @title Skeletal muscle differential analysis results
#'
#' @description Differential analysis of transcripts, proteins, phosphorylation
#'   sites, acetylation sites, redox sites, and metabolites measured in rat
#'   skeletal muscle (gastrocnemius and vastus lateralis).
#'
#' @usage
#' # Gastrocnemius
#' TRNSCRPT_GN_DA # Transcriptomics
#' PROT_GN_DA     # Proteomics
#' PHOSPHO_GN_DA  # Phosphoproteomics
#' ACETYL_GN_DA   # Acetyl proteomics
#' REDOX_GN_DA    # Redox proteomics
#' METAB_GN_DA    # Metabolomics
#'
#' # Vastus lateralis
#' TRNSCRPT_VL_DA # Transcriptomics
#' METAB_VL_DA    # Metabolomics
#'
#' @format A named list of 2 \code{data.frame} objects:
#'
#'   \describe{
#'     \item{"trained_vs_SED"}{All trained timepoints compared to their
#'     sex-matched sedentary control group. Example: F_1W - F_SED. Total of 8
#'     contrasts.}
#'     \item{"MvF"}{Males compared to females at SED and 8W timepoints.}
#'   }
#'
#'   In addition to all columns present in the \code{Biobase::fData} tables,
#'   each \code{data.frame} contains the following variables, most of which
#'   are added by \code{\link[limma:topTable]{limma::topTable}}:
#'
#'   \describe{
#'     \item{\code{contrast}}{factor; the contrast(s) being tested.}
#'     \item{\code{featureName}}{character; uniquely defines each feature that
#'     was tested.}
#'     \item{\code{logFC}}{numeric; difference in the weighted mean log\eqn{_2}
#'     values of the groups specified in the \code{contrast} column.}
#'     \item{\code{AveExpr}}{numeric; mean log\eqn{_2} relative abundance of
#'     each feature.}
#'     \item{\code{t}}{numeric; \acronym{LIMMA} moderated t-statistic.}
#'     \item{\code{z}}{numeric; Z-score.}
#'     \item{\code{P.Value}}{numeric; two-tailed p-value.}
#'     \item{\code{adj.P.Val}}{numeric; Benjamini-Hochberg adjusted p-value.
#'     P-values are adjusted across all trained vs. SED or male vs. female
#'     contrasts.}
#'     \item{\code{df.total}}{numeric; sum of residual and prior degrees of
#'     freedom.}
#'     \item{\code{signif_protein}}{logical; (PTM results only) whether the
#'     corresponding protein was differentially expressed.}
#'   }
#'
#' @details Code to reproduce the analyses is provided in the
#'   "differential_analysis" vignette.
#'
#' @seealso \code{\link[limma]{topTable}}
#'
#' @examples
#' # Number of differential features (FDR = 0.05)
#' f <- function(x) {
#'   lapply(x, function(.x) with(.x, table(contrast, adj.P.Val < 0.05)))
#' }
#'
#' ## Gastrocnemius
#' f(TRNSCRPT_GN_DA) # Transcripts
#' f(PROT_GN_DA)     # Proteins
#' f(PHOSPHO_GN_DA)  # Phosphosites
#' f(ACETYL_GN_DA)   # Acetylation sites
#' f(REDOX_GN_DA)    # Redox sites
#' f(METAB_GN_DA)    # Metabolites
#'
#' ## Vastus lateralis
#' f(TRNSCRPT_VL_DA) # Transcripts
#' f(METAB_VL_DA)    # Metabolites
#'
#' @keywords datasets
#'
#' @name SKM_DA

#' @rdname SKM_DA
#' @format NULL
#' @usage NULL
"TRNSCRPT_GN_DA"

#' @rdname SKM_DA
#' @format NULL
#' @usage NULL
"PROT_GN_DA"

#' @rdname SKM_DA
#' @format NULL
#' @usage NULL
"PHOSPHO_GN_DA"

#' @rdname SKM_DA
#' @format NULL
#' @usage NULL
"PHOSPHO_GN_NORM_DA"

#' @rdname SKM_DA
#' @format NULL
#' @usage NULL
"ACETYL_GN_DA"

#' @rdname SKM_DA
#' @format NULL
#' @usage NULL
"ACETYL_GN_NORM_DA"

#' @rdname SKM_DA
#' @format NULL
#' @usage NULL
"REDOX_GN_DA"

#' @rdname SKM_DA
#' @format NULL
#' @usage NULL
"REDOX_GN_NORM_DA"

#' @rdname SKM_DA
#' @format NULL
#' @usage NULL
"METAB_GN_DA"

#' @rdname SKM_DA
#' @format NULL
#' @usage NULL
"TRNSCRPT_VL_DA"

#' @rdname SKM_DA
#' @format NULL
#' @usage NULL
"METAB_VL_DA"


## GST results =================================================================

#' @title Molecular Signature Analysis Results
#'
#' @description Fast gene set enrichment analysis (FGSEA) and pre-ranked
#'   Correlation Adjusted MEan RAnk gene set testing (CAMERA) of gastrocnemius
#'   and vastus lateralis data. The gene sets available for testing, as well as
#'   dataset-specific gene sets, are detailed in \code{\link{rat_gene_sets}}.
#'   For metabolomics, RefMet chemical subclass sets were constructed from the
#'   "sub_class" column in the differential analysis results. Each -omics
#'   dataset was collapsed to the gene or RefMet metabolite/lipid ID level
#'   (metabolomics/lipidomics only) by selecting the most extreme z-score per
#'   contrast. See the vignettes for details.
#'
#' @usage
#'
#' ## FGSEA ----
#' # Gastrocnemius
#' TRNSCRPT_GN_FGSEA # Transcriptomics
#' PROT_GN_FGSEA     # Proteomics
#' ACETYL_GN_FGSEA   # Acetyl proteomics
#' REDOX_GN_FGSEA    # Redox proteomics
#' METAB_GN_FGSEA    # Metabolomics/lipidomics
#'
#' # Vastus lateralis
#' TRNSCRPT_VL_FGSEA # Transcriptomics
#' METAB_VL_FGSEA    # Metabolomics/lipidomics
#'
#' ## CAMERA-PR ----
#' # Gastrocnemius
#' TRNSCRPT_GN_CAMERA # Transcriptomics
#' PROT_GN_CAMERA     # Proteomics
#' ACETYL_GN_CAMERA   # Acetyl proteomics
#' REDOX_GN_CAMERA    # Redox proteomics
#' METAB_GN_CAMERA    # Metabolomics/lipidomics
#'
#' # Vastus lateralis
#' TRNSCRPT_VL_CAMERA # Transcriptomics
#' METAB_VL_CAMERA    # Metabolomics/lipidomics
#'
#'
#' @format Named lists of 2 \code{data.frame} objects:
#'
#' \describe{
#'   \item{"trained_vs_SED"}{All trained timepoints compared to their
#'   sex-matched sedentary control group. Example: F_1W - F_SED. Total of 8
#'   contrasts.}
#'   \item{"MvF"}{Males compared to females at SED and 8W timepoints.}
#' }
#'
#'   Each \code{data.frame} contains the following columns:
#'
#' \describe{
#'   \item{\code{"Contrast"}}{the comparison being made. Indicates which
#'   subset of the differential analysis results were used for testing.}
#'   \item{\code{"GeneSetID"}}{a unique identifier for each gene set (metabolomics and phosphoproteomics excluded). Used for cross-referencing and shortening descriptions for plotting.}
#'   \item{\code{"GeneSet"}}{a description of the gene set that was tested (metabolomics and phosphoproteomics excluded).}
#'   \item{\code{"GeneSetShort"}}{a shortened version of \code{"GeneSet"} used for plotting (metabolomics and phosphoproteomics only). If the description in \code{"GeneSet"} was above a certain length, the string is broken at the space between words and the unique \code{"GeneSetID"} is appended at the end.}
#'   \item{\code{"Subclass"}}{Metabolomics/lipidomics results only. A
#'   description of the RefMet chemical subclass set that was tested.}
#'   \item{\code{"Kinase"}}{Phosphoproteomics results only. The kinase that was
#'   tested.}
#'   \item{\code{"NGenesRat"}}{number of rat genes in the DA results for that
#'   contrast that match genes in \code{"GeneSet"}.}
#'   \item{\code{"NMetabolites"}}{Metabolomics/lipidomics results only. The
#'   number of RefMet metabolite/lipid identifiers from the DA results for that
#'   contrast that match a chemical subclass in \code{"Subclass"}.}
#'   \item{\code{"NSites"}}{Phosphoproteomics results only. The number of human
#'   sites from the DA results for that contrast that match sites in
#'   \code{"Kinase"}.}
#'   \item{\code{"NGenesHuman"}}{the size of the human gene set in MSigDB,
#'   before any filtering.}
#'   \item{\code{"NRatio"}}{\code{NGenesRat} / \code{NGenesHuman}. A measure of
#'   confidence that the set being tested is what is described by
#'   \code{"GeneSet"}.}
#'   \item{\code{"Direction"}}{direction of change of the gene set. Factor with
#'   levels "Up", "Down".}
#'   \item{\code{"NES"}}{FGSEA results only. The normalized enrichment score (NES). The true enrichment score (ES) divided by the mean of the permutation enrichment scores with the same sign.}
#'   \item{\code{"ZScore"}}{Z-score.}
#'   \item{\code{"PValue"}}{two-sided p-values.}
#'   \item{\code{"PAdj"}}{the Benjamini–Hochberg-adjusted p-value. P-values are
#'   adjusted across all gene sets and contrasts.}
#'   \item{\code{"LeadingEdge"}}{FGSEA only. A list of genes that contributed to
#'   the enrichment score. Usually, this is a subset of all members of the set.}
#' }
#'
#' @seealso \code{\link[fgsea]{fgsea}}, \code{\link[TMSig]{cameraPR.matrix}},
#'   \code{\link{molecular_signatures}}
#'
#' @examples
#' # Number of differential molecular signatures (FDR = 0.05)
#' f <- function(x) {
#'   lapply(x, function(.x) with(.x, table(Contrast, PAdj < 0.05)))
#' }
#'
#' ## FGSEA ----
#' # Gastrocnemius
#' #f1(TRNSCRPT_GN_FGSEA) # Transcriptomics
#' #f1(PROT_GN_FGSEA)     # Proteomics
#' #f1(PHOSPHO_GN_FGSEA)  # Phosphoproteomics
#' #f1(ACETYL_GN_FGSEA)   # Acetyl proteomics
#' #f1(REDOX_GN_FGSEA)    # Redox proteomics
#' #f1(METAB_GN_FGSEA)    # Metabolomics/lipidomics
#'
#' # Vastus lateralis
#' #f1(TRNSCRPT_VL_FGSEA) # Transcriptomics
#' #f1(METAB_VL_FGSEA)    # Metabolomics/lipidomics
#'
#' ## CAMERA-PR ----
#' # Gastrocnemius
#' #f1(TRNSCRPT_GN_CAMERA) # Transcriptomics
#' #f1(PROT_GN_CAMERA)     # Proteomics
#' #f1(PHOSPHO_GN_CAMERA)  # Phosphoproteomics
#' #f1(ACETYL_GN_CAMERA)   # Acetyl proteomics
#' #f1(REDOX_GN_CAMERA)    # Redox proteomics
#' #f1(METAB_GN_CAMERA)    # Metabolomics/lipidomics
#'
#' # Vastus lateralis
#' #f1(TRNSCRPT_VL_CAMERA) # Transcriptomics
#' #f1(METAB_VL_CAMERA)    # Metabolomics/lipidomics
#'
#'
#' @keywords datasets
#'
#' @name SKM_GST

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"TRNSCRPT_GN_FGSEA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"TRNSCRPT_GN_CAMERA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"PROT_GN_FGSEA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"PROT_GN_CAMERA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"ACETYL_GN_FGSEA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"ACETYL_GN_CAMERA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"ACETYL_GN_NORM_FGSEA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"ACETYL_GN_NORM_CAMERA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"PHOSPHO_GN_FGSEA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"PHOSPHO_GN_CAMERA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"PHOSPHO_GN_NORM_FGSEA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"PHOSPHO_GN_NORM_CAMERA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"REDOX_GN_FGSEA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"REDOX_GN_CAMERA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"REDOX_GN_NORM_FGSEA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"REDOX_GN_NORM_CAMERA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"METAB_GN_FGSEA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"METAB_GN_CAMERA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"TRNSCRPT_VL_FGSEA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"TRNSCRPT_VL_CAMERA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"METAB_VL_FGSEA"

#' @rdname SKM_GST
#' @format NULL
#' @usage NULL
"METAB_VL_CAMERA"


## Miscellaneous ===============================================================

#' @title Molecular Signatures
#'
#' @description Gene sets from the C5 and C8 collections of the Molecular
#'   Signatures Database (MSigDB) and phosphorylation sites grouped by their
#'   known kinases from PhosphoSitePlus. Used to summarize results of
#'   differential analysis performed on -omics datasets.
#'
#' @usage
#' human_gene_sets
#' rat_gene_sets
#' human_kinase_sets
#'
#' ## Dataset-specific sets ----
#' # Gastrocnemius
#' trnscrpt_GN_gene_sets  # Transcriptomics
#' prot_GN_gene_sets      # Proteomics
#' phospho_GN_kinase_sets # Phosphoproteomics
#' acetyl_GN_gene_sets    # Acetyl proteomics
#' redox_GN_gene_sets     # Redox proteomics
#' metab_GN_refmet_sets   # Metabolomics
#'
#' # Vastus lateralis
#' trnscrpt_VL_gene_sets  # Transcriptomics
#' metab_GN_gene_sets     # Metabolomics
#'
#'
#' @format Named lists of length 10472 (\code{human_gene_sets}), 9293
#'   (\code{rat_gene_sets}), or 319 (\code{human_kinase_sets}).
#'
#' @details Human gene sets were obtained from version 2023.2.Hs of the
#'   Molecular Signatures Database (MSigDB): specifically, those gene sets from
#'   the C5:GO subcollection (Gene Ontology terms) and the
#'   "RUBENSTEIN_SKELETAL_MUSCLE" sets from the C8 collection (Ashburner, 2000;
#'   Liberzon 2011, 2015; Rubenstein, 2020; The Gene Ontology Consortium, 2023).
#'
#'   Human genes were mapped to their rat orthologs using
#'   \code{\link[babelgene]{orthologs}} (v22.9; \code{\link{RAT_TO_HUMAN}}).
#'   Only those rat gene sets with at least 5 genes were retained.
#'
#'   Kinase–Substrate relationship data was obtained from
#'   "Kinase_Substrate_Dataset.gz", downloaded from PhosphoSitePlus
#'   (\acronym{PSP}) v6.7.1.1 (last modified Fri Nov 17 08:50:20 EST 2023):
#'   \url{https://www.phosphosite.org/staticDownloads.action}. Only human
#'   kinases and substrates were retained, since the data is substantially more
#'   comprehensive for humans than any other organism. Furthermore, only those
#'   kinases with at least 3 known substrates were kept. This results in
#'   phosphorylation sites grouped into 319 kinases.
#'
#'   Metabolomics RefMet subclass information was downloaded from
#'   \url{https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmetF_form.php}
#'   (accessed 2024-23-03).
#'
#'   Dataset-specific molecular signatures were generated by restricting each
#'   set to only those features measured in a particular tissue/ome combination
#'   with \code{\link[TMSig]{filterSets}} and only keeping those sets with at
#'   least 5 genes or metabolites/lipids or 3 human phosphorylation sites.
#'
#' @note PhosphoSitePlus data is not for commercial use.
#'
#' @references Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler,
#'   H., Cherry, J. M., Davis, A. P., Dolinski, K., Dwight, S. S., Eppig, J. T.,
#'   Harris, M. A., Hill, D. P., Issel-Tarver, L., Kasarskis, A., Lewis, S.,
#'   Matese, J. C., Richardson, J. E., Ringwald, M., Rubin, G. M., & Sherlock,
#'   G. (2000). Gene Ontology: Tool for the unification of biology. \emph{Nature
#'   Genetics, 25}(1), 25–29. \url{https://doi.org/10.1038/75556}
#'
#'   Liberzon, A., Subramanian, A., Pinchback, R., Thorvaldsdóttir, H., Tamayo,
#'   P., & Mesirov, J. P. (2011). Molecular signatures database (MSigDB) 3.0.
#'   \emph{Bioinformatics, 27}(12), 1739–1740.
#'   \url{https://doi.org/10.1093/bioinformatics/btr260}
#'
#'   Hornbeck, P. V., Zhang, B., Murray, B., Kornhauser, J. M., Latham, V., &
#'   Skrzypek, E. (2015). PhosphoSitePlus, 2014: mutations, PTMs and
#'   recalibrations. \emph{Nucleic acids research, 43}, D512–D520.
#'   \url{https://doi.org/10.1093/nar/gku1267}
#'
#'   Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P.,
#'   & Tamayo, P. (2015). The Molecular Signatures Database (MSigDB) hallmark
#'   gene set collection. \emph{Cell systems, 1}(6), 417–425.
#'   \url{https://doi.org/10.1016/j.cels.2015.12.004}
#'
#'   Rubenstein, A. B., Smith, G. R., Raue, U., Begue, G., Minchev, K.,
#'   Ruf-Zamojski, F., Nair, V. D., Wang, X., Zhou, L., Zaslavsky, E., Trappe,
#'   T. A., Trappe, S., & Sealfon, S. C. (2020). Single-cell transcriptional
#'   profiles in human skeletal muscle. \emph{Scientific reports, 10}(1), 229.
#'   \url{https://doi.org/10.1038/s41598-019-57110-6}
#'
#'   Martin, F. J., Amode, M. R., Aneja, A., Austine-Orimoloye, O., Azov, A. G.,
#'   Barnes, I., Becker, A., Bennett, R., Berry, A., Bhai, J., Bhurji, S. K.,
#'   Bignell, A., Boddu, S., Branco Lins, P. R., Brooks, L., Ramaraju, S. B.,
#'   Charkhchi, M., Cockburn, A., Da Rin Fiorretto, L., … Flicek, P. (2023).
#'   Ensembl 2023. \emph{Nucleic Acids Research, 51}(D1), D933–D941.
#'   \url{https://doi.org/10.1093/nar/gkac958}
#'
#'   The Gene Ontology Consortium, Aleksander, S. A., Balhoff, J., Carbon, S.,
#'   Cherry, J. M., Drabkin, H. J., Ebert, D., Feuermann, M., Gaudet, P.,
#'   Harris, N. L., Hill, D. P., Lee, R., Mi, H., Moxon, S., Mungall, C. J.,
#'   Muruganugan, A., Mushayahama, T., Sternberg, P. W., Thomas, P. D., …
#'   Westerfield, M. (2023). The Gene Ontology knowledgebase in 2023.
#'   \emph{Genetics, 224}(1), iyad031.
#'   \url{https://doi.org/10.1093/genetics/iyad031}
#'
#' @keywords datasets
#'
#' @name molecular_signatures

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"human_gene_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"rat_gene_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"SET_TO_ID"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"acetyl_GN_NORM_gene_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"acetyl_GN_gene_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"human_kinase_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"metab_GN_refmet_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"metab_VL_refmet_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"phospho_GN_NORM_kinase_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"phospho_GN_kinase_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"prot_GN_gene_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"redox_GN_NORM_gene_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"redox_GN_gene_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"trnscrpt_GN_gene_sets"

#' @rdname molecular_signatures
#' @format NULL
#' @usage NULL
"trnscrpt_VL_gene_sets"

#' @title Rat-to-human gene conversion table
#'
#' @description Output of \code{\link[babelgene]{orthologs}} (babelgene v22.9).
#'
#' @examples
#' head(RAT_TO_HUMAN)
#'
#' @keywords datasets
"RAT_TO_HUMAN"

#' @title Rat-to-human gene conversion table for phophosites
#' @description Output of \code{\link[babelgene]{orthologs}} (babelgene v22.9).
#' @keywords datasets
"RAT_TO_HUMAN_SITE"


#' @title Redundant gene sets
#'
#' @description Gene sets that were removed following use of
#'   \code{\link[TMSig]{clusterSets}} with default parameters. These tend to be
#'   subsets of those sets that were retained for analysis.
#'
#' @usage
#' TRNSCRPT_GN_REDUNDANT_SETS
#' PROT_GN_REDUNDANT_SETS
#' PHOSPHO_GN_REDUNDANT_SETS
#' ACETYL_GN_REDUNDANT_SETS
#' REDOX_GN_REDUNDANT_SETS
#'
#' TRNSCRPT_VL_REDUNDANT_SETS
#'
#' @format an object of class \code{data.frame} with 5 columns:
#'
#' \describe{
#'   \item{"GeneSetID"}{character; a unique gene set ID. Mainly used for fast
#'   cross-referencing with the FGSEA or CAMERA-PR results.}
#'   \item{"GeneSet"}{character; the gene set that was retained for the
#'   analysis. Usually, the largest set in the cluster. This column does not
#'   contain all gene sets tested—only those that appeared in a cluster with
#'   other sets.}
#'   \item{"SimilarGeneSets"}{character; the gene set that was excluded for analysis.}
#'   \item{"JaccardCoef"}{numeric; the Jaccard similarity coefficient. A value
#'   of 1 indicates that "GeneSet" and "SimilarGeneSets" are aliased, at least
#'   after they had been restricted to only those genes measured in the
#'   dataset.}
#'   \item{"OverlapCoef"}{numeric; the overlap similarity coefficient. A value
#'   of 1 indicates that "SimilarGeneSets" is a subset of "GeneSet", or they are
#'   aliased, at least after they had been restricted to only those genes
#'   measured in the dataset.}
#' }
#'
#' @keywords datasets
#'
#' @name SKM_REDUNDANT_SETS

#' @rdname SKM_REDUNDANT_SETS
#' @format NULL
#' @usage NULL
"TRNSCRPT_GN_REDUNDANT_SETS"

#' @rdname SKM_REDUNDANT_SETS
#' @format NULL
#' @usage NULL
"PROT_GN_REDUNDANT_SETS"

#' @rdname SKM_REDUNDANT_SETS
#' @format NULL
#' @usage NULL
"ACETYL_GN_REDUNDANT_SETS"

#' @rdname SKM_REDUNDANT_SETS
#' @format NULL
#' @usage NULL
"REDOX_GN_REDUNDANT_SETS"

#' @rdname SKM_REDUNDANT_SETS
#' @format NULL
#' @usage NULL
"TRNSCRPT_VL_REDUNDANT_SETS"

#' @rdname SKM_REDUNDANT_SETS
#' @format NULL
#' @usage NULL
"ACETYL_GN_NORM_REDUNDANT_SETS"

#' @rdname SKM_REDUNDANT_SETS
#' @format NULL
#' @usage NULL
"REDOX_GN_NORM_REDUNDANT_SETS"



#' @title Differential analysis incidence tables
#'
#' @usage
#' INCIDENCE_DA
#'
#' @format A named list of two \code{data.frame} objects:
#'
#' \describe{
#'   \item{"trained_vs_SED"}{All trained timepoints compared to their
#'   sex-matched sedentary control group. Example: F_1W - F_SED. Total of 8
#'   contrasts.}
#'   \item{"MvF"}{Males compared to females at SED and 8W timepoints.}
#' }
#'
#' Each \code{data.frame} contains the following columns:
#'
#' \describe{
#'   \item{tissue}{the type of skeletal muscle. Either "GN" (gastrocnemius) or
#'   "VL" (vastus lateralis).}
#'   \item{ome}{factor; the type of -omics dataset. Either "TRNSCRPT", "PROT",
#'   "ACETYL", "REDOX", or "METAB".}
#'   \item{featureName}{character; the feature ID. May be a protein, transcript
#'   ID, PTM site, or metabolite/lipid.}
#'   \item{gene_symbol}{character; the corresponding HGNC gene symbol. Not
#'   applicable for metabolites.}
#'   \item{key_female}{character; ("trained_vs_SED" results only) a unique key
#'   indicating the female trained vs. SED contrasts in which the feature was
#'   differentially abundant.}
#'   \item{key_male}{character; ("trained_vs_SED" results only) a unique key
#'   indicating the male trained vs. SED contrasts in which the feature was
#'   differentially abundant.}
#'   \item{key}{character; ("MvF" results only) a unique key indicating the
#'   contrasts in which the feature was differentially abundant.}
#' }
#'
#' The remaining columns are for the contrasts (8 for trained vs. SED, 2 for
#' male vs. female). A value of 1 indicates that gene set or feature
#' significantly differed (BH-adjusted p-value < 0.05) between groups, while a
#' value of 0 indicates they did not. These values are collapsed to form the
#' \code{key_female} and \code{key_male} columns for more convenient filtering.
#'
#' @keywords datasets
#'
#' @name INCIDENCE_DA
