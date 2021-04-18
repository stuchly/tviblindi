#' Compute filtration and generate boundary matrices
#'
#' This method adds new slots \code{filtration}, \code{boundary} and \code{reduced_boundary} to an object of class \code{tviblindi}.
#' In filtration, a simplicial complex is built over the \code{data} point cloud and a filtration value is associated with each simplex.
#'
#' @param tv \code{tviblindi}-class object with \code{kNN}
#' @param method string: type of simplicial complex. Either \code{witness} or \code{alpha}. Default value is \code{witness}. Computing \code{alpha} complex for dimensionality > 5 will probably fail (memory requirements scale badly with dimension of input data)
#' @param k integer: nearest neighbour count for finding nearest vertices in \code{denoised} (or \code{data} if \code{denoised} is not available) to SOM cluster centers in \code{codes}
#' @param alpha number: relaxation parameter for witness complex. Default value is \code{NULL}, whereby the mean distance between matched vertices in negihbourhoods of \code{k} is taken
#' 
#' @references
#' \insertRef{Edelsbrunner2008}{tviblindi}
#' 
#' \insertRef{Bates2019}{tviblindi}
#' 
#' \insertRef{Rouvreau2020}{tviblindi}
#' 
#' \insertRef{Maria2020}{tviblindi}
#' 
#' \insertRef{Kachanovich2020}{tviblindi}
#' 
#' \insertRef{Bauer2014}{tviblindi}
#' 
#' \insertRef{Fasey2019}{tviblindi}
#' 
#' \insertRef{Eddelbuettel2013}{tviblindi}
#'
#' @export
Filter.tviblindi <- function(
  tv,
  method = 'witness',
  k = 30,
  alpha = NULL
) {
  if (!is.character(method) || (!method[1] %in% c('witness', 'alpha')))
    stop(paste0('Invalid filtration method: ', method))
  if (is.null(tv$codes))
    stop('Data not clustered')
  
  use_original <- is.null(tv$denoised)
  
  if (method == 'witness') {
    .msg('Constructing witness complex')
    matched_vertices <- FNN::get.knnx(
      data = tv$codes,
      query = if (use_original) { tv$data } else { tv$denoised },
      k = k
    )
    if (is.null(alpha)) {
      alpha <- mean(matched_vertices$nn.dist[, k])
      .msg_alt(paste0('alpha = ', alpha))
    }
    index_list <- split(matched_vertices$nn.index, seq(nrow(matched_vertices$nn.index)))
    distance_list <- split(matched_vertices$nn.dist,  seq(nrow(matched_vertices$nn.index)))
    tv$filtration <- witness_from_distances_cliques(index_list, distance_list, alpha2 = alpha, maxdimension = 1)
    gc(verbose = FALSE)
    tv$filtration <- create_k_skeleton(coordinates = tv$codes, filtration = tv$filtration, k = 2)
  } else if (method == 'alpha') {
    .msg('Constructing alpha complex')
    tv$filtration <- TDA::alphaComplexFiltration(tv$codes, printProgress = TRUE)
  }
  gc(verbose = FALSE)
  
  tv$boundary <- boundary_matrix(tv$filtration)
  tv$reduced_boundary <- reduced_boundary_matrix(tv$boundary)
  
  tv$filtration_method <- paste0(
    method,
    ' (', if (method == 'witness') { paste0('alpha = ', alpha, ', ') } else { '' }, paste0('k = ', k), ')'
  )
  
  invisible(tv)
}

Filter <- function(tv, method, k, alpha) UseMethod('Filter', tv)

#' Create k-skeleton of a simplicial complex
#'
#' Provided a simplicial complex filtration and corresponding point coordinates, this function  generates a \code{k}-skeleton of the complex.
#' 
#' @param filtration list: filtration of a simplicial complex. Must contain definition (vertex indices) of each simplex by indices in \code{cmplx}. \code{values} slot (with filtration values) is disregarded, as filtration values must are re-computed from scratch
#' @param coordinates numeric matrix of coordinates: rows of matrix correspond to vertices, number of columns is number of dimensions
#' @param k integer: skeleton \code{k} (dimensionality)
#' @param type character: type of filtration. Only allowed value now is '\code{alpha}'
#' 
#' @return list; new filtration, with \code{cmplx} (simplex indices) and \code{values} (filtration values for each simplex)
#' 
#' @export
create_k_skeleton <- function(
  filtration,
  coordinates,
  k,
  type = 'alpha'
) {
  which_cmplx <- which(lapply(filtration$cmplx, length) <= (k + 1))
  cmplx <- filtration$cmplx[which_cmplx]
  d <- unlist(lapply(cmplx, length)) - 1 # dimensions
  
  vals <- vector(mode = 'numeric')
  if (type == 'alpha') {
    vals <- alpha_complex_filtration_values(cmplx, coordinates)
  } else {
    stop('Invalid value of "type"')
  }
  ord <- order(vals, d)
  
  result <- list()
  
  result$cmplx <- cmplx[ord]
  result$values <- vals[ord]  
  
  return(result)
}

#' Compute alpha-complex filtration values
#'
#' Provided a simplicial complex filtration and corresponding point coordinates, this function calculates the corresponding filtration values for each simplex.
#' \code{filtration} must be a list containing the slot \code{cmplx}, which is a list of integer vectors: 1-based vertex indices.
#' These define simplices in the alpha-complex filtration.
#' 
#' @param filtration list: filtration of a simplicial complex. Must contain definition (vertex indices) of each simplex by indices in \code{cmplx}. \code{values} slot (with filtration values) is disregarded, as filtration values must are re-computed from scratch
#' @param coordinates numeric matrix of coordinates: rows of matrix correspond to vertices, number of columns is number of dimensions
#' 
#' @return list: new filtration, with \code{cmplx} (simplex indices) and \code{values} (filtration values for each simplex)
#'
#' @export
alpha_complex_filtration_values <- function(
  filtration,
  coordinates
) {
  filtration.u <- sort(unique(unlist(filtration)))
  return(alpha_complex_filtration_values_C(filtration, coordinates, filtration.u))
}

#' Compute boundaries of a filtration
#'
#' Provided a simplicial complex filtration, \code{boundary_matrix} returns a list of boundary simplices (\code{cmplx}) and unmodified filtration values (\code{values}).
#' \code{filtration} must be a list containing definition (vertex indices) of each simplex by indices in \code{cmplx} and filtration values of each simplex in \code{values}.
#' 
#' @param filtration list: filtration of a simplicial complex. Must contain description of each simplex by indices in \code{cmplx} and associated filtration values in \code{values}
#' 
#' @return list with slots \code{cmplx} (boundary simplex indices) and \code{values} (filtration values for each simplex)
#'
#' @export
boundary_matrix <- function(
  filtration
) {
  b <- build_boundary_matrix(filtration$cmplx)
  list(
    cmplx = b,
    values = filtration$values
  )
}

#' Reduce boundary matrix of a filtration
#'
#' Provided a boundary matrix corresponding to a simplicial complex filtration, \code{reduced_boundary_matrix} returns a list with slots \code{boundary} (mapping of simplices to their boundaries), \code{nonzero_col}, \code{low} and \code{dim} (see reference below for description of the algorithm).
#' 
#' @param boundary list containing the slot \code{cmplx}, which is a list of vectors (1-based simplex indices), and the slot \code{values} containing filtration values of each simplex
#' 
#' @return list with slots \code{boundary}, \code{nonzero_col}, \code{low} and \code{dim}.
#'
#' @references
#' \insertRef{Edelsbrunner2008}{tviblindi}
#'
#' @export
reduced_boundary_matrix <- function(
  boundary
) {
  values <- boundary$values
  boundary$cmplx <- lapply(boundary$cmplx, function(x) if (length(x)==0) return(integer(0)) else return(x-1))
  boundary <- phat_boundary(boundary, 1)
  gc(verbose = FALSE)
  boundary$values <- values
  boundary
}
