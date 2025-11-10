#' @title Robust Wrapper for Linear Modeling with LIMMA
#'
#' @description A \acronym{LIMMA} wrapper for -omics data in the form of
#'   \code{\link[Biobase]{eSet}} or similar objects, including a separate
#'   pipeline for RNA-Seq count data.
#'
#' @param object object inheriting from class \code{\link[Biobase]{eSet}}.
#' @param model.str character; model formula of the form \code{"~ a + b"} or
#'   \code{"~ 0 + a + b"} (no intercept), where \code{"a"}, and \code{"b"} are
#'   predictors (usually columns of \code{pData(object)}.
#' @param contrasts character; (optional) one or more contrasts of the form
#'   \code{"varA - varB"}, where \code{var} is a predictor included in
#'   \code{model.str}. Contrasts involving interactions are not allowed. See
#'   \code{\link{paircomp}} for assistance with creating these contrasts.
#' @param trend logical, character; should an intensity-dependent trend be
#'   allowed for the prior variance? If \code{FALSE}, then the prior variance is
#'   constant. Alternatively, \code{trend} can be the name of a column in
#'   \code{fData(object)}, which will be used as the covariate for the prior
#'   variance. Passed to \code{\link[limma]{eBayes}}.
#' @param robust logical; should the estimation of \code{df.prior} and
#'   \code{var.prior} be made robust to hyper- or hypo-variable features? Passed
#'   to \code{\link[limma]{eBayes}}.
#' @param var.group character; name of a column in \code{pData(object)} that
#'   groups samples for \code{\link[limma]{arrayWeights}}. If \code{NULL}
#'   (default), samples will be weighted equally.
#' @param block character; name of a column in \code{pData(object)} specifying a
#'   blocking factor used in \code{\link[limma]{duplicateCorrelation}}. If
#'   \code{NULL} (default), \code{duplicateCorrelation} is not used.
#' @param plot logical; whether to show model summary/diagnostic plots.
#'
#' @section Model Specification:
#'
#'   If at least one predictor in \code{model.str} is categorical, exclusion of
#'   the intercept merely changes the parameterization of the model from an
#'   \strong{effects} model to a \strong{means} model.
#'
#'   \strong{Effects} model:
#'
#'   \deqn{y_{ij} = \mu_i + \epsilon_{ij}}
#'
#'   \strong{Means} model:
#'
#'   \deqn{y_{ij} = \mu + \tau_i + \epsilon_{ij}}
#'
#'   where \eqn{y_{ij}} is the \eqn{ij}th observation, \eqn{\mu_i} is the mean
#'   of the \eqn{i}th level of the categorical predictor, \eqn{\mu} is the
#'   overall mean, \eqn{\tau_i} is the \eqn{i}th treatment effect (difference
#'   between the mean of the \eqn{i}th level and the overall mean), and
#'   \eqn{\epsilon_{ij}} is the random error.
#'
#'   A no-intercept model is required if \code{contrasts} are provided.
#'
#'   If \code{model.str} is a no-intercept model with all numeric predictors, an
#'   intercept will be recommended with a warning, since this makes the (often
#'   incorrect) assumption that the regression line should pass through the
#'   origin.
#'
#' @section General Details:
#'
#'   If \code{trend} is the name of a column in \code{fData(object)}, missing
#'   values are not allowed.
#'
#'   If there are fewer than 10 features, a warning will be issued, both
#'   \code{robust} and \code{trend} will be set to \code{FALSE}, and samples
#'   will be weighted equally.
#'
#' @section Modeling (RNA-Seq) Count Data:
#'
#'   If \code{exprs(object)} is a matrix of counts, the voom pipeline for
#'   RNA-Seq count data will be used. In that case, \code{trend} will be set to
#'   \code{FALSE}, since \code{\link[edgeR]{voomLmFit}} incorporates the
#'   mean-variance trend into the precision weights. Low-count features will be
#'   removed (with a warning) by \code{\link[edgeR]{filterByExpr}}.
#'
#'   If there are fewer than 10 features, a warning will be issued, \code{trend}
#'   will be set to \code{FALSE}, and samples will be weighted equally.
#'
#' @return
#'
#' An object of class \code{\link[limma:MArrayLM-class]{MArrayLM}}.
#'
#' @references
#'
#' Anders, S, Huber, W (2010). Differential expression analysis for sequence
#' count data. \emph{Genome Biology} 11, R106.
#'
#' Bullard JH, Purdom E, Hansen KD, Dudoit S. (2010) Evaluation of statistical
#' methods for normalization and differential expression in mRNA-Seq
#' experiments. \emph{BMC Bioinformatics} 11, 94.
#'
#' Chen Y, Lun ATL, and Smyth, GK (2016). From reads to genes to pathways:
#' differential expression analysis of RNA-Seq experiments using Rsubread and
#' the edgeR quasi-likelihood pipeline. \emph{F1000Research} 5, 1438.
#' \url{https://f1000research.com/articles/5-1438}
#'
#' Law, CW, Chen, Y, Shi, W, Smyth, GK (2014). Voom: precision weights unlock
#' linear model analysis tools for RNA-seq read counts. \emph{Genome Biology}
#' 15, R29. See also the Preprint Version at
#' \url{http://www.statsci.org/smyth/pubs/VoomPreprint.pdf} incorporating some
#' notational corrections.
#'
#' Lun, ATL, Smyth, GK (2017). No counts, no variance: allowing for loss of
#' degrees of freedom when assessing biological variability from RNA-seq data.
#' Statistical Applications in Genetics and Molecular Biology 16(2), 83-93.
#' \href{https://www.doi.org/10.1515/sagmb-2017-0010}{doi:10.1515/sagmb-2017-0010}
#'
#' Liu, R, Holik, AZ, Su, S, Jansz, N, Chen, K, Leong, HS, Blewitt, ME,
#' Asselin-Labat, M-L, Smyth, GK, Ritchie, ME (2015). Why weight? Modelling
#' sample and observational level variability improves power in RNA-seq
#' analyses. \emph{Nucleic Acids Research 43}, e97.
#' \href{https://www.doi.org/10.1093/nar/gkv412}{doi:10.1093/nar/gkv412}
#'
#' Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016). Robust
#' hyperparameter estimation protects against hypervariable genes and improves
#' power to detect differential expression. \emph{Annals of Applied Statistics}
#' 10, 946-963. \url{http://projecteuclid.org/euclid.aoas/1469199900}
#'
#' Ritchie, M. E., Diyagama, D., Neilson, van Laar, R., J., Dobrovic, A.,
#' Holloway, A., and Smyth, G. K. (2006). Empirical array quality weights in the
#' analysis of microarray data. \emph{BMC Bioinformatics} \strong{7}, 261.
#' \href{https://doi.org/10.1186/1471-2105-7-261}{doi:10.1186/1471-2105-7-261}
#'
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, and Smyth GK (2015). limma
#' powers differential expression analyses for RNA-sequencing and microarray
#' studies. \emph{Nucleic Acids Research} 43, e47.
#' \url{http://nar.oxfordjournals.org/content/43/7/e47}
#'
#' Robinson MD, Oshlack A (2010). A scaling normalization method for
#' differential expression analysis of RNA-seq data. \emph{Genome Biology} 11,
#' R25.
#'
#' Smyth, G. K. (2002). An efficient algorithm for REML in heteroscedastic
#' regression. \emph{Journal of Computational and Graphical Statistics}
#' \strong{11}, 836-847. \url{https://gksmyth.github.io/pubs/remlalgo.pdf}
#'
#' Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing
#' differential expression in microarray experiments. \emph{Statistical
#' Applications in Genetics and Molecular Biology}, \strong{3}, No. 1, Article
#' 3. https://gksmyth.github.io/pubs/ebayes.pdf
#'
#' Smyth, G. K., Michaud, J., and Scott, H. (2005). The use of within-array
#' replicate spots for assessing differential expression in microarray
#' experiments. \emph{Bioinformatics 21}(9), 2067-2075.
#' \url{http://bioinformatics.oxfordjournals.org/content/21/9/2067}
#' Preprint with corrections: \url{https://gksmyth.github.io/pubs/dupcor.pdf}
#'
#' @seealso
#'
#' \code{\link[limma]{limmaUsersGuide}}, \code{\link[limma]{lmFit}},
#' \code{\link{paircomp}}
#'
#' @examples
#' library(Biobase)
#'
#' # Simulate ExpressionSet object
#' set.seed(0)
#' nsamples <- 12L
#' nfeatures <- 1000L
#' mat_dimnames <- list(paste0("feature.", seq_len(nfeatures)),
#'                      paste0("sample.", seq_len(nsamples)))
#' mat <- matrix(rnorm(nsamples * nfeatures), ncol = nsamples,
#'               dimnames = mat_dimnames)
#'
#' p_data <- data.frame(group = rep(LETTERS[1:2], each = nsamples / 2),
#'                      covariate = rpois(n = nsamples, lambda = 3),
#'                      row.names = colnames(mat))
#' p_data <- AnnotatedDataFrame(p_data)
#'
#' obj <- ExpressionSet(assayData = mat, phenoData = p_data)
#'
#' # Intercept model, categorical predictor
#' fit1 <- limmaFit(obj, model.str = "~ group")
#' head(fit1$coefficients)
#'
#' # No-intercept model, categorical predictor
#' fit2 <- limmaFit(obj, model.str = "~ 0 + group")
#' head(fit2$coefficients)
#'
#' # Intercept model, categorical and numeric predictors
#' fit3 <- limmaFit(obj, model.str = "~ group + covariate")
#' head(fit3$coefficients)
#'
#' # No-intercept model, categorical and numeric predictors
#' fit4 <- limmaFit(obj, model.str = "~ 0 + group + covariate")
#' head(fit4$coefficients)
#'
#' # Include contrasts (requires no-intercept model)
#' fit5 <- limmaFit(obj, model.str = "~ 0 + group",
#'                  contrasts = "groupA - groupB")
#' head(fit5$coefficients) # coefficients same as fit5$coefficients[, "groupB"]
#'
#' @importFrom limma is.fullrank nonEstimable arrayWeights duplicateCorrelation
#'   lmFit makeContrasts contrasts.fit eBayes plotMDS plotSA
#' @importFrom edgeR DGEList filterByExpr normLibSizes voomLmFit
#' @importFrom Biobase exprs pData fData fData<- featureNames
#' @importFrom stats terms model.matrix
#' @importFrom graphics abline barplot par
#'
#' @export limmaFit


if (!require("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")

library(limma)
library(edgeR)
library(Biobase)
library(stats)
library(graphics)

limmaFit <- function(object,
                     model.str,
                     contrasts,
                     trend = TRUE,
                     robust = TRUE,
                     var.group = NULL,
                     block = NULL,
                     plot = FALSE)
{
  ## Check arguments -----------------------------------------------------------
  stopifnot(inherits(object, what = "eSet"))
  stopifnot(is.logical(robust))
  .checkArgs(object, model.str, trend, var.group, block)

  plot <- plot & interactive()

  # Add "featureName" column to fData. This will be included in the output of
  # limmaDEA.
  fData(object)[["featureName"]] <- featureNames(object)
  # Reorder columns
  fData_cols <- c("featureName",
                  setdiff(colnames(fData(object)), "featureName"))
  fData(object) <- fData(object)[, fData_cols, drop = FALSE]

  # Remove rows and columns with all NA values
  object <- .dropNAVectors(object)

  # Convert model.str to formula for model.matrix
  model.formula <- eval(parse(text = model.str), envir = pData(object))

  # Extract model information
  model.terms <- terms(model.formula)
  # all.vars excludes higher-order terms and functions of variables
  model.vars <- all.vars(model.formula, functions = FALSE, unique = TRUE)
  term.labels <- attr(model.terms, "term.labels")

  # Check that all lower-order terms are columns of pData
  missing_terms <- setdiff(model.vars, colnames(pData(object)))
  if (length(missing_terms) > 0L)
    stop("The following `model.str` terms are not columns of pData(object): ",
         paste(missing_terms, collapse = ", "))

  # If contrasts are provided, check that they are valid
  if (!missing(contrasts))
    .validateContrasts(object, contrasts, model.terms)

  # Suggest an intercept if all model.str terms are numeric
  all_terms_numeric <- all(
    sapply(pData(object)[, model.vars, drop = FALSE],
           function(pred_i) is.numeric(pred_i))
  )

  no_intercept <- attr(model.terms, which = "intercept") == 0

  if (no_intercept & all_terms_numeric)
    warning("An intercept is recommended when all terms in ",
            "`model.str` are numeric.")

  # Suggest an intercept if an interaction is present
  if (no_intercept & any(grepl(":", term.labels)))
    warning("An intercept is recommended when interaction terms are ",
            "present in `model.str`.")

  # Remove samples with one or more NA model terms
  NA_terms <- apply(is.na(pData(object)[, model.vars, drop = FALSE]), 1, any)

  if (all(NA_terms))
    stop("All samples were removed due to one or more NA model terms.")

  if (any(NA_terms)) {
    warning("Removing the following samples by index due to one or more NA ",
            "model terms: ", paste(which(NA_terms), collapse = ", "),
            ". These samples will not contribute to the analysis.")
    object <- object[, !NA_terms]
  }

  ### Create design matrix
  design <- model.matrix(model.formula)
  .checkFullRankDesign(design, model.formula)

  ### Linear modeling ----------------------------------------------------------

  # Automatically detect count data
  count_data <- isTRUE(all.equal(exprs(object), floor(exprs(object))))
  if (count_data) {
    ## RNA-seq (count) data pipeline -------------------------------------------
    message("Using the voom pipeline for RNA-Seq count data.")

    # Convert to DGEList
    dge <- DGEList(counts = exprs(object),
                   samples = pData(object),
                   genes = fData(object))

    # Filter to transcripts with sufficiently large counts
    keep <- filterByExpr(dge, design = design, group = NULL)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    num_features <- sum(keep) # number of features remaining

    if (num_features == 0L)
      stop("All features were removed by edgeR::filterByExpr.")

    if (any(!keep))
      message(sum(!keep),
              " low count features were removed by edgeR::filterByExpr(). ",
              num_features, " features remaining.")

    # Calculate library size normalization factors
    dge <- normLibSizes(dge, method = "TMM")

    # If there are fewer than 10 features, do not weight samples or make the
    # global variance robust to the presence of outliers.
    if (num_features < 10L) {
      if (!is.null(var.group)) {
        warning("Fewer than 10 features. Samples will be weighted equally.")
        var.group <- NULL
      }
      if (robust) {
        warning("Fewer than 10 features. eBayes() will use robust=FALSE.")
        robust <- FALSE
      }
    } else if (!is.null(var.group)) {
      var.group <- dge$samples[[var.group]]
    }

    # Uses arrayWeights and duplicateCorrelation if var.group and block are not
    # NULL, respectively.
    fit <- voomLmFit(counts = dge, design = design,
                     block = block, var.group = var.group,
                     plot = plot, save.plot = plot, keep.EList = TRUE)

    if (!isFALSE(trend)) {
      warning("voom incorporates the mean-variance trend into the ",
              "precision weights, so eBayes() will use trend=FALSE.")
      trend <- FALSE
    }

    if (plot) {
      # Set weights for barplot later
      if (is.null(fit$targets$sample.weight)) {
        weights <- rep(1, ncol(dge))
      } else {
        weights <- fit$targets$sample.weight
      }
      names(weights) <- seq_along(weights)

      # MDS plot (log2 normalized counts per million)
      # TODO allow custom labels? Defaults to sample names.
      plotMDS(fit$EList$E)
    }

  } else {
    ## Non RNA-seq data pipeline -----------------------------------------------

    # MDS plot
    # TODO allow custom labels? Defaults to sample names.
    if (plot)
      plotMDS(exprs(object))

    # Default sample weights (samples are weighted equally)
    weights <- rep(1, ncol(object))
    names(weights) <- seq_along(weights)

    num_features <- nrow(object)

    if (num_features < 10) {
      warning("Fewer than 10 features. Samples will be weighted equally ",
              "and eBayes() will use trend=FALSE, robust=FALSE.")
      trend <- robust <- FALSE
    } else if (!is.null(var.group)) {
      weights <- arrayWeights(object, design, method = "genebygene",
                              var.group = pData(object)[[var.group]])
    }

    if (!is.logical(trend))
      trend <- fData(object)[[trend]]

    # Estimate correlation within blocks (necessary for repeated measures,
    # where the blocks are biological replicates)
    if (!is.null(block)) {
      # message("Estimating the intra-block correlation...")
      block <- eval(parse(text = block), envir = pData(object))
      dupecor <- duplicateCorrelation(object, design, block = block,
                                      weights = weights)$consensus.correlation

      fit <- lmFit(object, design, weights = weights,
                   block = block, correlation = dupecor)
    } else {
      fit <- lmFit(object, design, weights = weights)
    }
  }

  # Contrasts
  if (!missing(contrasts)) {
    contrast.matrix <- makeContrasts(contrasts = contrasts,
                                     levels = design)
    fit <- contrasts.fit(fit, contrasts = contrast.matrix)
  }

  ## Empirical Bayes moderation
  fit.smooth <- eBayes(fit, trend = trend, robust = robust)

  # Diagnostic plots
  if (plot) {
    oldpar <- par(mfrow = c(1, 2))
    on.exit(par(oldpar))

    # Sample weights
    barplot(weights, space = 0, xlab = "Sample Index", ylab = "Weight")
    abline(a = 1, b = 0, lty = "longdash", col = "red")

    # Mean-variance trend (should be a flat line if voom was used)
    plotSA(fit.smooth)
  }

  return(fit.smooth)
}



## Helper functions ------------------------------------------------------------

# Check limmaFit arguments.

#' @importFrom Biobase pData fData
.checkArgs <- function(object,
                       model.str,
                       trend = TRUE,
                       var.group = NULL,
                       block = NULL)
{
  # Check that model.str is a string
  if (!is.character(model.str) | length(model.str) > 1L)
    stop("`model.str` must be a length 1 character string.")

  if (!grepl("^[ ]*~", model.str))
    stop("`model.str` should begin with a \"~\" followed by relevant ",
         "model terms. Unexpected `model.str` format.")

  # Check that var.group is a column in pData, if not NULL
  if (!is.null(var.group)) {
    if (is.null(pData(object)[[var.group]]))
      stop("`var.group` must be NULL or the name of a column in pData(object).")
  }

  # Check that block is a column in pData, if not NULL
  if (!is.null(block)) {
    if (is.null(pData(object)[[block]]))
      stop("`block` must be NULL or the name of a column in pData(object).")
  }

  # Check that trend is a column in pData with no NA values, if not logical
  if (!is.logical(trend)) {
    if (is.null(fData(object)[[trend]])) {
      stop("If `trend` is not logical, it must be the name of a column",
           "in fData(object).")
    } else if (any(is.na(fData(object)[[trend]]))) {
      stop("NA values not allowed in fData(object)[[trend]].")
    }
  }
}


# Remove rows and columns with all missing expression values.

#' @importFrom Biobase exprs
.dropNAVectors <- function(object) {
  x <- is.na(exprs(object))

  if (all(x))
    stop("No non-missing values in exprs(object).")

  keep <- lapply(1:2, function(i) {
    dim_i <- c(" rows", " columns")[i]
    idx_i <- apply(x, i, all) # TRUE = all NA
    n <- sum(idx_i)

    if (n)
      warning("There are ", n, dim_i, " with ",
              "all missing values. These will be removed.")

    return(!idx_i) # TRUE = not all NA
  })

  object <- object[keep[[1]], keep[[2]]]

  return(object)
}


# If contrasts are provided, check the model specification.

#' @importFrom Biobase pData
.validateContrasts <- function(object,
                               contrasts,
                               model.terms)
{
  if (!is.character(contrasts))
    stop("`contrasts` must be a vector of one or more character strings ",
         "specifying comparisons between levels of a term in `model.str`.")

  term.labels <- attr(model.terms, "term.labels")

  # If contrasts are provided, a no-intercept model is required to ensure they
  # are specified correctly
  if (attr(model.terms, which = "intercept") == 1) {
    model.str <- paste("~ 0 +", paste(term.labels, collapse = " + "))

    stop("If contrasts are provided, please specify a no-intercept model, ",
         "like so: model.str = ", dQuote(model.str, FALSE))
  }
}


# Ensure design matrix is full-rank. If any columns of the design matrix are
# linearly dependent, output those columns in an error message.

#' @importFrom limma is.fullrank nonEstimable
#' @importFrom stats terms
.checkFullRankDesign <- function(design,
                                 model.formula) {
  if (!is.fullrank(design)) {
    # Map nonestimable design matrix columns to the model.str terms
    design_idx <- which(colnames(design) %in% nonEstimable(design))
    term_num <- attr(design, "assign")[design_idx]
    nonest_terms <- unique(attr(terms(model.formula), "term.labels")[term_num])
    # Ignore interactions
    nonest_terms <- nonest_terms[!grepl(":", nonest_terms)]

    if (length(nonest_terms) > 0L)
      stop("The design matrix is not full rank. The following model terms ",
           "are linearly dependent on previous terms: ",
           paste(nonest_terms, collapse = ", "))
  }
}
