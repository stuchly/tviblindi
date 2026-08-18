#' Dimensionality reduction via sparse diffusion maps
#'
#' Built-in "diffuse" \code{DimRed} method. Registered automatically at package load; see
#' \code{\link{register_dimred_method}} for the method-plugin contract this follows.
#'
#' @param x tviblindi class object; must already have \code{x$KNN} computed (see \code{\link{KNN}}).
#' @param data numeric matrix; the resolved input data (raw or denoised). Unused directly - diffuse
#' maps operate on \code{x$KNN} instead - but present for signature consistency with other methods.
#' @param neigen integer (default 2); number of eigenvectors to compute.
#' @param t double (default 0); time parameter; if \code{t==0} multi-time scale is used (geometric sum).
#' @param ... unused; absorbs any extra arguments DimRed's dispatch may forward.
#'
#' @return numeric matrix of reduced coordinates.
#'
#' @export
dimred_diffuse <- function(x, data, neigen = 2, t = 0, ...) {
  if (is.null(x$KNN)) stop("Compute KNN graph first.")
  sparse.diffuse(
    sparse.Laplacian.construct(knn.raw2adj(x$KNN)),
    neigen = neigen,
    t = t
  )$X
}
