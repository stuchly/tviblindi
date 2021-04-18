.rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x / as.vector(sm)
}

.upsample <- function(
  data,
  knn_idcs,
  idcs,
  n,
  concentration = 0.1,
  labels = NULL
) {
  lab <- out <- NULL
  for (i in idcs) {
    e <- matrix(.rdirichlet(n, rep(concentration, ncol(knn_idcs))), byrow = TRUE, ncol = n)
    out <- rbind(out, t(t(data[knn_idcs[i, ], ]) %*% e))
    if (!is.null(labels))
      lab <- c(lab, rep(labels[i], n))
  }
  list(out = out, lab = lab)
}

#' Create a new \code{tviblindi} object with upsampled data
#'
#' \code{Upsample} creates a new \code{tviblindi} object, taking the associated input expression matrix with its labels and sampling new additional points.
#' Dirichlet distribution over neighbours of each point is used to sample the new points.
#' Any slots besides \code{data} and \code{labels} will be cleared by this.
#'
#' @param tv \code{tviblindi}-class object with a k-nearest-neighbour graph computed
#' @param n integer: how many points shall be sampled for each upsampled point in original data. Default value is 5
#' @param k integer: number of nearest neighbours of each original point to use for fitting the Dirichlet distribution to draw new point from. Default value is 30
#' @param idx_labels integer or string: index of the vector of labels to use for re-labelling the data. Default value is 1
#' @param idcs integer or string vector: row indices of points in \code{tv$data} to which upsampling shall be applied, or a vector of labels of these points. Default value is \code{NULL} (this reverts to all points in \code{tv$data})
#' @param concentration float: alpha parameter of the Dirichlet distribution. Default value is 0.1
#'
#' @export
Upsample <- function(
  tv,
  n = 5,
  k = 30,
  idx_labels = 1,
  idcs = NULL,
  concentration = 0.1
) {
  if (is.null(idcs))
    idcs <- 1:nrow(tv$data)
  if (is.character(idcs))
    idcs <- which(tv$labels[[idx_labels]] %in% unique(idcs))
  k <- min(k, ncol(tv$kNN$idcs))
  
  res <- .upsample(tv$data, tv$kNN$idcs[, 1:k], idcs, n, concentration, as.character(tv$labels[[idx_labels]]))
  labels <- c(as.character(tv$labels[[idx_labels]]), res$lab)
  data <- rbind(tv$data, res$out)
  
  tv <- tviblindi(data = data, labels = labels)
  tv
}

