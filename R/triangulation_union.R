# triangulation (filtration) union

unifyFiltrations <- function(
  filt, # list of filtrations, must include 'cmplx' and 'values'
  sortInput=TRUE
) {
  cmplx <- unlist(lapply(filt, function(x) x$cmplx), recursive=FALSE) # merge cmplx entries
  values <- unlist(lapply(filt, function(x) x$values), recursive=FALSE) # merge cmplx entries

  inds <- unique_simplex_inds(cmplx, sortInput) # in src/triangulation_union.cpp

  cmplx <- cmplx[inds]
  values <- values[inds]

  d <- unlist(lapply(cmplx,length)) # dimensions (+1)
  ord <- order(values, d) # order by filtration values, then by dimension

  # produce merged filtration:
  jointFilt <- list("cmplx"=cmplx[ord], "values"=values[ord], "increasing"=TRUE)

  warning("filtration values might be degenerated, use create_k_skeleton to recompute them")

  return(jointFilt)
}
