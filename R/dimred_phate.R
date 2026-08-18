#' Dimensionality reduction via PHATE
#'
#' Built-in "phate" \code{DimRed} method. Registered automatically at package load; see
#' \code{\link{register_dimred_method}} for the method-plugin contract this follows.
#'
#' @param x tviblindi class object; unused, present for signature consistency with other methods.
#' @param data numeric matrix; the resolved input data (raw or denoised).
#' @param ndim integer (default 2); dimension of reduced data.
#' @param ... additional arguments forwarded to \code{phateR::phate} (e.g. \code{knn}, \code{decay},
#' \code{t}, \code{gamma}, \code{npca}).
#'
#' @return numeric matrix of reduced coordinates.
#'
#' @export
dimred_phate <- function(x, data, ndim = 2, ...) {
  if (!requireNamespace("phateR", quietly = TRUE)) stop("install package 'phateR' first")
  layout <- as.matrix(phateR::phate(data, ndim = ndim, ...))
  rownames(layout) <- NULL
  layout
}
