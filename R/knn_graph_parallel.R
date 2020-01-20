knn.adj.raw.parallel <- function(X, K, metric="euclid") {
  d <- 0
  if (metric == "cosine")
    d <- 1
  return(openmp_knn_C(X, K, d))
}
