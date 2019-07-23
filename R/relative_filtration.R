# relative filtration

relativeFiltration <- function(indices, filtration, sortOutput=FALSE) {
  if ((length(indices) == 0) || (length(filtration) == 0)) {
    message("indices must be a non-empty integer vector, filtration a non-empty list of integer vectors") } else {
    relative_filtration(indices, filtration, sortOutput) }
}
