#' @title Median-MAD-normalize expression data
#'
#' @description Normalize a matrix of log\eqn{_2}-transformed expression data.
#'   Median-center columns and optionally scale by the mean absolute deviations
#'   (MADs).
#'
#' @param x a matrix with features (genes, proteins, etc.) as rows and samples
#'   as columns.
#' @param mad logical; whether to divide the median-centered samples by the mean
#'   absolute deviations (MADs) and then multiply by the mean of the MADS.
#'   Default \code{FALSE} only median-centers the samples.
#'
#' @returns output of \code{\link[base]{scale}}. The centered, optionally scaled
#'   matrix.
#'
#' @examples
#' set.seed(0)
#' x <- matrix(rnorm(200, sd = 2), nrow = 20, ncol = 10)
#' boxplot(x) # before normalization
#'
#' x1 <- normalize_data(x, mad = FALSE)
#' boxplot(x1) # median-center samples only
#'
#' x2 <- normalize_data(x, mad = TRUE)
#' boxplot(x2) # median-MAD normalize
#'
#' # Sample medians and MADs
#' attr(x2, "scaled:center") # medians
#' attr(x2, "scaled:scale") # mads
#'
#' @importFrom stats median mad
#'
#' @export normalize_data

normalize_data <- function(x, mad = FALSE) {
  # Sample medians
  medians <- apply(x, 2, median, na.rm = TRUE)

  if (mad) {
    # Sample MADs
    mads <- apply(x, 2, mad, na.rm = TRUE)
    s <- mean(mads)

    x <- scale(x, center = medians, scale = mads)
    x <- x * s
  } else {
    # Subtract sample medians only
    x <- scale(x, center = medians, scale = FALSE)
  }

  return(x)
}
