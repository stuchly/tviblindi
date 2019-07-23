#' Create k-skeleton of a simplicial complex
#'
#' Provided a simplicial complex filtration and corresponding point coordinates, \code{create_k_skeleton} generates a skeleton of the given complex in \code{k} dimensions.
#' 
#' @param filtration list; filtration of a simplicial complex. Must contain description of each simplex by indices in \code{cmplx}. \code{values} slot is disregarded, as filtration values must be computed from scratch.
#' @param coordinates matrix of coordinates. Rows of matrix correspond to vertices, columns correspond to dimensions.
#' @param k skeleton dimensionality.
#' @param type type of filtration: 'alpha' or 'rips'. Default is 'alpha'.
#' @param threshold Rips threshold, used it type is 'rips'.
#' 
#' @return list; new filtration, with \code{cmplx} (simplex indices) and \code{values} (filtration values for each simplex).
#' 
#'
#' @export
create_k_skeleton <- function(filtration, coordinates, k, type="alpha", threshold=0) # filtration values for a k-skeleton
{
  which_cmplx   <- which(lapply(filtration$cmplx, length)<=(k+1))
  cmplx         <- filtration$cmplx[which_cmplx]
  d             <- unlist(lapply(cmplx,length))-1 # dimensions
  
  vals          <- vector(mode="numeric")
  if (type=="alpha") {
    vals <- alpha_complex_filtration_values(cmplx, coordinates)
  } else if (type=="rips") {
    vals <- rips_complex_filtration_values(cmplx, coordinates, threshold)
  } else {
    warning("create_k_skeleton param type must be valid, see documentation")
    return(-1)
  }
  ord           <- order(vals, d)
  
  result        <- list()
  
  result$cmplx  <- cmplx[ord]
  result$values <- vals[ord]  

  return(result)
}
