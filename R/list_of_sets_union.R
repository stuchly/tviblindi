# union of lists of integer vectors

unifyListsOfSets <- function(
  L, # list of lists of integer vectors
  sortInput=TRUE
) {
  sets <- unlist(L, recursive=FALSE) # merge cmplx entries
  
  inds <- unique_simplex_inds(sets, sortInput) # in src/triangulation_union.cpp
  
  return(sets[inds])
}
