tiny_tviblindi_with_knn <- function() {
  set.seed(2)
  data <- matrix(rnorm(60 * 4), ncol = 4)
  labels <- sample(c("a", "b"), 60, replace = TRUE)
  x <- tviblindi(data = data, labels = list(pop = labels))
  KNN(x, K = 10)
  x
}

test_that("dimred_diffuse requires a KNN graph", {
  x <- tviblindi(data = matrix(rnorm(20 * 3), ncol = 3), labels = list(pop = rep("a", 20)))
  expect_error(DimRed(x, method = "diffuse"), "Compute KNN graph first")
})

test_that("DimRed(method='diffuse') matches calling the diffusion-map primitives directly", {
  # Regression check for the refactor: the pre-refactor DimRed.tviblindi computed the
  # diffuse layout inline as this exact expression. dimred_diffuse() should be a pure
  # extraction, so results must be bit-for-bit identical for the same x$KNN and RNG state
  # (ARPACK's Lanczos start vector is randomized, so eigenvector sign - not magnitude -
  # varies run to run unless the seed is pinned immediately before each call).
  x <- tiny_tviblindi_with_knn()

  set.seed(123)
  DimRed(x, method = "diffuse", neigen = 3, t = 0)
  got <- x$layout[[1]]

  set.seed(123)
  expected <- sparse.diffuse(
    sparse.Laplacian.construct(knn.raw2adj(x$KNN)),
    neigen = 3,
    t = 0
  )$X

  expect_equal(got, expected)
  expect_equal(dim(got), c(nrow(x$data), 3))
})
