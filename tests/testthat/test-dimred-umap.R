tiny_tviblindi <- function() {
  set.seed(3)
  data <- matrix(rnorm(50 * 4), ncol = 4)
  labels <- sample(c("a", "b"), 50, replace = TRUE)
  tviblindi(data = data, labels = list(pop = labels))
}

test_that("DimRed(method='umap') returns a finite [N,2] layout", {
  skip_if_not_installed("uwot")
  x <- tiny_tviblindi()
  set.seed(42)
  DimRed(x, method = "umap", n_neighbors = 5)

  got <- x$layout[[1]]
  expect_equal(dim(got), c(nrow(x$data), 2))
  expect_true(all(is.finite(got)))
  expect_null(rownames(got)) # rownames explicitly stripped, matching pre-refactor behavior
})

test_that("extra ... arguments are forwarded through to uwot::umap", {
  skip_if_not_installed("uwot")
  x <- tiny_tviblindi()
  expect_error(DimRed(x, method = "umap", not_a_real_uwot_arg = TRUE))
})
