#' Reduce noise in high-dimensional biological data in a \code{tviblindi} object
#'
#' This method adds a denoised expression matrix to a \code{tviblindi} object.
#' Each vertex' coordinates are replaced by the average of its nearest neighbours' coordinates.
#' This is repeated \code{n_iter} times.
#' This modification hopefully brings the data points closer to the original manifold and removes technical noise.
#'
#' @param tv \code{tviblindi}-class object
#' @param k integer: number of neighbours to use for averaging the coordinates. Default value is \code{30}
#' @param n_iter integer: number of iterations. Default value is \code{1}
#'
#' @export
Denoise.tviblindi <- function(
  tv,
  k = 30,
  n_iter = 1
) {
  if (is.null(tv$kNN))
    stop('k-NN graph needs to be produced first')
  
  max_k <- ncol(tv$kNN$idcs)
  
  if (k > max_k) {
    k <- max_k
    warning('k is larger than nearest neighbour count in kNN matrix')
  }
  
  idcs <- tv$kNN$idcs[, 2:k]
  
  tv$denoised <- denoise_matrix(tv$data, idcs, k - 1, n_iter)
  colnames(tv$denoised) <- colnames(tv$data)
  tv$denoising_iterations <- tv$denoising_iterations + n_iter
  
  .msg(paste0('Data denoised in ', n_iter, ' iteration', if (n_iter > 1) 's' else ''))
  
  gc(verbose = FALSE)
  
  invisible(tv)
}

Denoise <- function(tv, k, n_iter) UseMethod('Denoise', tv)

#' Set cell-of-origin
#'
#' This method adds a new cell-of-origin to a \code{tviblindi} object.
#' A point with the chosen index is set as cell-of-origin, in which simulated random walks start.
#' Alternatively, a point in the data is found such that it is closest to the centroid of a chosen annotated population and picked as a cell-of-origin.
#'
#' @param tv \code{tviblindi}-class object
#' @param label integer or string: either label of population of origin or index of the cell-of-origin
#' @param name string: name for the particular cell-of-origin. Default value is '\code{default}'
#'
#' @export
SetOrigin.tviblindi <- function(
  tv,
  label,
  name = 'default'
) {
  if (is.null(tv$data))
    stop('No expression data contained in the tviblindi analysis object')
  if (length(label) != 1)
    stop('A single character label or integer index must be given as label')
  
  if (is.integer(label)) {
    if (label < 1 || label > nrow(tv$data)) stop('Origin index out of bounds')
    
    if (is.null(tv$origin)) {
      tv$origin <- list()
      tv$origin_label <- list()
    }
    if (name %in% names(tv$origin))
      warning(paste0('Origin "', name, '" already exists (index = ', tv$origin[[name]], '). Overwriting'))
    
    tv$origin[[name]] <- label
    tv$origin_label[[name]] <- as.character(label)
    
  } else {
    stopifnot(is.character(label))
    stems <- which(tv$labels[[name]] == label)
    if (length(stems) == 0)
      stop(paste0('No events with label ', label, ' were found'))
    
    if (is.null(tv$origin)) {
      tv$origin <- list()
      tv$origin_label <- list()
    }
    
    if (name %in% names(tv$origin[[name]]))
      warning(paste0('Origin "', name, '" already exists (index = ', tv$origin[[name]], '). Overwriting'))
    
    tv$origin[[name]] <-
      stems[which.min(rowSums(t(t(tv$data[stems, , drop = FALSE]) - colMeans(tv$data[stems, , drop = FALSE]))^2))]
    tv$origin_label[[name]] <- as.character(label)
  }
  
  .msg(paste0('Origin set to vertex ', tv$origin[[name]]))
  invisible(tv)
}

SetOrigin <- function(tv, label, name) UseMethod('SetOrigin', tv)

#' Compute a k-nearest-neighbours matrix
#'
#' This method adds a new slot \code{kNN} to an object of class \code{tviblindi}.
#' For exact k-NN, a multi-threaded version of ball-tree nearest neighbour search is used.
#' For the approximate version, the \code{annoy} algorithm is used.
#' 
#' A list \code{kNN} with slots \code{idcs} and \code{dist}, containing matrices with nearest neighbour indices and distances (respectively), with each row corresponding to a row of \code{data}, is produced.
#'
#' @param tv \code{tviblindi}-class object
#' @param k integer: nearest neighbour count
#' @param approximate logical: whether to apply an approximate algorithm (faster). Default value is \code{TRUE}
#' @param n_trees integer: number of look-up trees to be constructed by \code{annoy}, if applicable. Default value is \code{150}
#'
#' @references
#' \insertRef{Nl2008}{tviblindi}
#' 
#' \insertRef{Ulyanov2016}{tviblindi}
#' 
#' \insertRef{Bernhardsson2020}{tviblindi}
#'
#' @export
ConstructkNNG.tviblindi <- function(
  tv,
  k           = 100,
  approximate = TRUE,
  n_trees     = 150
) {
  stopifnot(!is.null(tv))
  if (is.null(tv$data))
    stop('No expression data contained in the tviblindi analysis object')
  tv$kNN <- if (approximate) kNN_annoy(tv$data, k, trees = n_trees) else kNN_exact(tv$data, k)
  
  gc(verbose = FALSE)
  invisible(tv)
}

ConstructkNNG <- function(tv, k, approximate, n_trees) UseMethod('ConstructkNNG', tv)

#' Construct a k-NN matrix using the 'annoy' algorithm
#'
#' \code{kNN_annoy} takes a coordinate matrix \code{data}, a nearest neighbour count \code{k} and a look-up tree count \code{trees}.
#' It finds the \code{k} nearest neighbours to each vertex defined in \code{data} using an approximate algorithm.
#'
#' @param data numeric matrix: coordinates (one data point per each row)
#' @param k integer: nearest neighbour count
#' @param trees integer: number of look-up trees to generate by \code{annoy}. The higher the number of trees, the higher the accuracy. Default value is \code{150}
#'
#' @return a list \code{kNN} with slots \code{idcs} and \code{dist} containing matrices with nearest neighbour indices and distances, respectively, with each row corresponding to a row in \code{data}
#'
#' @references
#' \insertRef{Nl2008}{tviblindi}
#' 
#' \insertRef{Ulyanov2016}{tviblindi}
#' 
#' \insertRef{Bernhardsson2020}{tviblindi}
#'
#' @export
kNN_annoy <- function(
  data,
  k,
  trees = 150
) {
  res <- knn_annoy(data, k, trees)
  list(
    idcs = do.call(rbind, res[[1]]),
    dist = do.call(rbind, res[[2]])
  )
}
