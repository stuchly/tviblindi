#' Compute k-means or SOM clusters for triangulation of data
#'
#' This method adds new slots \code{clusters} and \code{codes} to an object of class \code{tviblindi}.
#' These contain (1) the mapping of vertices in \code{data} onto nodes in a newly trained k-means model or self-organising map (SOM) and (2) cluster centers, respectively.
#'
#' @param tv \code{tviblindi}-class object with \code{kNN}
#' @param method string: clustering method specification, either '\code{som}' or '\code{kmeans}'. Default value is '\code{k-means}'
#' @param kmeans_algorithm string: algorithm for k-means clustering. One of '\code{Hartigan-Wong}', '\code{Lloyd}', '\code{Forgy}', '\code{MacQueen}'. Default value is '\code{Hartigan-Wong}'
#' @param k integer: cluster count for k-means clustering. Default value is 625
#' @param xdim integer: width of Kohonen map computed by \code{FlowSOM} (if \code{method == 'som'}). Default value is 25
#' @param ydim integer: height of Kohonen map computed by \code{FlowSOM} (if \code{method == 'som'}). Default value is 25
#'
#' @references
#' \insertRef{VanGassen2015}{tviblindi}
#'
#' @export
Cluster.tviblindi <- function(
  tv,
  method = 'kmeans',
  kmeans_algorithm = 'Hartigan-Wong',
  k = 625,
  xdim = 25,
  ydim = 25
) {
  use_denoised <- !is.null(tv$denoised)
  if (!use_denoised)
    .msg_alt_bad('Using non-denoised data')
  
  if (method == 'kmeans') {
    n_clusters <- k
    .msg(paste0('Sampling ', n_clusters, ' points as initial cluster centers'))
    pts <- do.call(rbind, internal_sample_points(if (use_denoised) tv$denoised else tv$data, n_clusters)$landmarks)
    gc(verbose = FALSE)
    
    res <- stats::kmeans(
      if (use_denoised) tv$denoised else tv$data,
      centers = pts, algorithm = kmeans_algorithm
    )
    tv$clusters       <- res$cluster
    tv$codes          <- res$centers
    tv$cluster_method <- method
    tv$k              <- k
    tv$som_xgrid      <- NULL
    tv$som_ygrid      <- NULL
  } else if (method == 'som') {
    n_clusters <- xdim * ydim
    .msg(paste0('Sampling ', n_clusters, ' points as initial cluster centers'))
    pts <- do.call(rbind, internal_sample_points(if (use_denoised) tv$denoised else tv$data, n_clusters)$landmarks)
    gc(verbose = FALSE)
    
    .msg(paste0('Constructing a ', xdim, ' by ', ydim, ' SOM grid'))
    som <- FlowSOM::SOM(
      if (use_denoised) tv$denoised else tv$data,
      codes = pts,
      xdim  = xdim,
      ydim  = ydim
    )
    tv$clusters       <- som$mapping[, 1]
    tv$codes          <- som$codes
    tv$cluster_method <- method
    tv$k              <- NULL
    tv$som_xgrid      <- xdim
    tv$som_ygrid      <- ydim
  } else {
    stop('Invalid clustering method')
  }
  invisible(tv)
}

Cluster <- function(tv, xdim, ydim) UseMethod('Cluster', tv)