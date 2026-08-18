tiny_tviblindi <- function() {
  set.seed(7)
  data <- matrix(rnorm(50 * 5), ncol = 5)
  labels <- sample(c("a", "b"), 50, replace = TRUE)
  tviblindi(data = data, labels = list(pop = labels))
}

# "randproj" is not a built-in method; it's registered here at runtime the same way a
# third-party user/script would, via the public register_dimred_method() API - see the
# dimred-methods skill / R/dimred_registry.R for the plugin contract.
randproj <- function(x, data, k = 2, ...) {
  d <- ncol(data)
  dirs <- matrix(rnorm(d * k), nrow = d, ncol = k)
  dirs <- apply(dirs, 2, function(v) v / sqrt(sum(v^2)))
  data %*% dirs
}

register_dimred_method("randproj", randproj, overwrite = TRUE)

test_that("DimRed(method='randproj') returns a finite [N,2] layout by default", {
  x <- tiny_tviblindi()
  set.seed(11)
  DimRed(x, method = "randproj")

  got <- x$layout[[1]]
  expect_equal(names(x$layout), "1_randproj")
  expect_equal(dim(got), c(nrow(x$data), 2L))
  expect_true(all(is.finite(got)))
  expect_true(is.matrix(got))
})

test_that("randproj's k= argument is forwarded through DimRed's ...", {
  x <- tiny_tviblindi()
  set.seed(12)
  DimRed(x, method = "randproj", k = 4)

  got <- x$layout[[1]]
  expect_equal(dim(got), c(nrow(x$data), 4L))
})

test_that("randproj projects onto unit vectors (columns of the projection are consistent with the data's own scale)", {
  x <- tiny_tviblindi()
  set.seed(13)
  DimRed(x, method = "randproj", k = 3)

  got <- x$layout[[1]]
  # each output column is a linear combination of the input row with unit-norm weights, so
  # its values must be bounded by the row's Euclidean norm (Cauchy-Schwarz)
  row_norms <- sqrt(rowSums(x$data^2))
  expect_true(all(abs(got) <= row_norms + 1e-8))
})
