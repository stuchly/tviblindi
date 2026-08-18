#' Dimensionality reduction via UMAP
#'
#' Built-in "umap" \code{DimRed} method. Registered automatically at package load; see
#' \code{\link{register_dimred_method}} for the method-plugin contract this follows.
#'
#' @param x tviblindi class object; unused, present for signature consistency with other methods.
#' @param data numeric matrix; the resolved input data (raw or denoised).
#' @param ... additional arguments forwarded to \code{uwot::umap}.
#'
#' @return numeric matrix of reduced coordinates.
#'
#' @export
dimred_umap <- function(x, data, ...) {
  if (!requireNamespace("uwot", quietly = TRUE)) stop("install package 'uwot' first")
  layout <- uwot::umap(data, verbose = TRUE, ...)
  rownames(layout) <- NULL
  layout
}
