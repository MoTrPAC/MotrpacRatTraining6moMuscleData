#' @title Robust Wrapper for Differential Analysis with LIMMA
#'
#' @description A wrapper function for \code{\link[limma]{topTable}}. Especially
#'   useful when performing moderated t-tests for more than one coefficient or
#'   contrast.
#'
#' @param fit an \code{\link[limma:MArrayLM-class]{MArrayLM}} object. For
#'   example, the output of \code{\link{limmaFit}}.
#' @param coef indices or names of columns of \code{fit$coefficients} to test.
#' @param test character; the type of statistical test to perform for each
#'   feature. Either \code{"t"} (moderated t-tests) or \code{"F"} (moderated
#'   F-tests).
#' @param adjust.method character; one of \code{p.adjust.methods}. Default is
#'   \code{"BH"}.
#' @param adjust.globally logical; whether to perform global p-value adjustment
#'   for sets of related contrasts.
#'
#' @returns Output of \code{\link[limma]{topTable}} with additional columns
#'   \code{"contrast"} and \code{"df.total"}.
#'
#' @details
#'
#' See \code{\link{limmaFit}} to easily create
#' \code{\link[limma:MArrayLM-class]{MArrayLM}} objects for \code{fit}.
#'
#' If \code{adjust.globally=TRUE}, the entire vector of p-values will be
#' adjusted together. Only do this if sets of related contrasts are being
#' tested. Otherwise, p-values are adjusted separately by contrast using
#' \code{adjust.method}. By default \code{adjust.method = "BH"}, which uses the
#' method of Benjamini & Hochberg (1995) to control the false discovery rate.
#'
#' @references
#'
#' Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate:
#' a practical and powerful approach to multiple testing. \emph{Journal of the
#' Royal Statistical Society Series B}, \strong{57}, 289–300.
#' doi:\href{https://doi.org/10.1111/j.2517-6161.1995.tb02031.x}{10.1111/j.2517-6161.1995.tb02031.x}.
#'
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, and Smyth GK (2015). limma
#' powers differential expression analyses for RNA-sequencing and microarray
#' studies. \emph{Nucleic Acids Research} 43, e47.
#' \url{http://nar.oxfordjournals.org/content/43/7/e47}
#'
#' @seealso
#'
#' \code{\link[limma]{limmaUsersGuide}}, \code{\link{limmaFit}},
#' \code{\link[limma]{topTable}}
#'
#' @importFrom data.table rbindlist
#' @importFrom limma topTable
#' @importFrom stats p.adjust p.adjust.methods
#'
#' @export limmaDEA

limmaDEA <- function(fit,
                     coef = NULL,
                     test = c("t", "F"),
                     adjust.method = "BH",
                     adjust.globally = FALSE)
{
  # Check input
  if (!"MArrayLM" %in% class(fit))
    stop("`fit` must be an object of class \"MArrayLM\".")

  test <- match.arg(test, choices = c("t", "F"))
  adjust.method <- match.arg(adjust.method, choices = p.adjust.methods)

  # Coefficient and contrast matrices
  coef.matrix <- fit$coefficients
  contrast.matrix <- fit[["contrasts"]]

  if (is.null(coef)) {
    if (is.null(contrast.matrix)) {
      test <- "F"
    } else {
      coef <- colnames(contrast.matrix)
    }
  } else if (is.numeric(coef)) {
    # convert indices to strings
    coef <- colnames(coef.matrix)[coef]
  }

  # Moderated F-tests
  if (test == "F") {
    if (length(coef) == 1L) {
      warning("test='F', but there is only one coefficient. Using test='t'.")
      test <- "t"
    } else {
      out <- topTable(fit,
                      coef = coef,
                      number = Inf,
                      adjust.method = adjust.method,
                      sort.by = "none")
      out$df.total <- fit$df.total

      # Rename coefficient columns
      coef_idx <- seq(ncol(fit$genes) + 1,
                      which(colnames(out) == "AveExpr") - 1)
      colnames(out)[coef_idx] <- paste0("coef.", colnames(out)[coef_idx])
    }
  }

  # Moderated t-tests
  if (test == "t") {
    out <- lapply(coef, function(coef_i) {
      x <- topTable(fit,
                    coef = coef_i,
                    number = Inf,
                    adjust.method = adjust.method,
                    sort.by = "none")
      x$df.total <- fit$df.total
      return(x)
    })
    names(out) <- coef
    out <- rbindlist(out, idcol = "contrast")
    out$contrast <- factor(out$contrast, levels = coef)

    if (adjust.globally)
      out$adj.P.Val <- p.adjust(out$P.Value, method = adjust.method)
  }

  rownames(out) <- NULL
  return(out)
}

