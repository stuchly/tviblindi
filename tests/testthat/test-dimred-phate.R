tiny_tviblindi <- function() {
  set.seed(3)
  data <- matrix(rnorm(50 * 4), ncol = 4)
  labels <- sample(c("a", "b"), 50, replace = TRUE)
  tviblindi(data = data, labels = list(pop = labels))
}

test_that("DimRed(method='phate') returns a finite [N,2] layout", {
  skip_if_not_installed("phateR")
  x <- tiny_tviblindi()
  set.seed(42)
  DimRed(x, method = "phate", knn = 5)

  got <- x$layout[[1]]
  expect_equal(dim(got), c(nrow(x$data), 2))
  expect_true(all(is.finite(got)))
  expect_null(rownames(got))
})

test_that("extra ... arguments are forwarded through to phateR::phate", {
  skip_if_not_installed("phateR")
  x <- tiny_tviblindi()
  expect_error(DimRed(x, method = "phate", not_a_real_phate_arg = TRUE))
})
