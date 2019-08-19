#' Generate boundary matrix object
#'
#' Provided a list of simplices, \code{build_boundaryC} makes list of preceding simplices which are boundaries. This is, effectively, a compact representation of a boundary matrix.
#'
#' @param dFilt list; filtration of a simplicial complex. Contains description of each simplex by indices in \code{cmplx} and filtration values in \code{values}.
#'
#' @return A list of \code{cmplx} (boundaries) and \code{values} (same as in \code{dFilt}).
#'
#'
#' @export
build_boundaryC<-function(dFilt){
  db<-build_boundary_C(dFilt$cmplx)
  return(list(cmplx=db,values=dFilt$values))
}
