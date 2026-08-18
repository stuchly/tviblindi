tiny_data <- function(n = 20, d = 4, seed = 9) {
  set.seed(seed)
  matrix(rnorm(n * d), ncol = d)
}

test_that("viz_data defaults to data when not supplied", {
  data <- tiny_data()
  x <- tviblindi(data = data, labels = list(pop = rep("a", nrow(data))))
  expect_identical(x$viz_data, data)
})

test_that("viz_data can be supplied explicitly at construction (dense)", {
  data <- tiny_data()
  viz <- tiny_data(seed = 10)
  x <- tviblindi(data = data, labels = list(pop = rep("a", nrow(data))), viz_data = viz)
  expect_identical(x$viz_data, viz)
})

test_that("viz_data can be supplied explicitly at construction (sparse)", {
  data <- tiny_data()
  viz <- Matrix::Matrix(tiny_data(seed = 11), sparse = TRUE)
  x <- tviblindi(data = data, labels = list(pop = rep("a", nrow(data))), viz_data = viz)
  expect_true(methods::is(x$viz_data, "Matrix"))
  expect_equal(dim(x$viz_data), dim(viz))
})

test_that("viz_data row-count mismatch errors at construction", {
  data <- tiny_data()
  viz <- tiny_data(n = 5)
  expect_error(
    tviblindi(data = data, labels = list(pop = rep("a", nrow(data))), viz_data = viz)
  )
})

test_that("Set_viz_data replaces viz_data in place and returns invisibly", {
  data <- tiny_data()
  x <- tviblindi(data = data, labels = list(pop = rep("a", nrow(data))))
  viz <- tiny_data(seed = 12)
  ret <- Set_viz_data(x, viz)
  expect_identical(x$viz_data, viz)
  expect_false(identical(x$viz_data, data))
  expect_true(is.environment(ret)) # invisible(x); still returns x
})

test_that("Set_viz_data row-count mismatch errors and leaves x unchanged", {
  data <- tiny_data()
  x <- tviblindi(data = data, labels = list(pop = rep("a", nrow(data))))
  before <- x$viz_data
  expect_error(Set_viz_data(x, tiny_data(n = 3)))
  expect_identical(x$viz_data, before)
})

test_that("Set_viz_data accepts a sparse matrix", {
  data <- tiny_data()
  x <- tviblindi(data = data, labels = list(pop = rep("a", nrow(data))))
  viz <- Matrix::Matrix(tiny_data(seed = 13), sparse = TRUE)
  Set_viz_data(x, viz)
  expect_true(methods::is(x$viz_data, "Matrix"))
})

test_that("DownSample.tviblindi resamples viz_data in lockstep with data", {
  data <- tiny_data(n = 40)
  x <- tviblindi(data = data, labels = list(pop = rep(c("a", "b"), 20)))
  viz <- tiny_data(n = 40, seed = 14)
  Set_viz_data(x, viz)
  KNN(x, K = 10)
  DownSample(x, N = 15, K = 8)
  expect_equal(nrow(x$viz_data), nrow(x$data))
  # row identity preserved: viz_data's surviving rows still match the pre-downsample viz rows
  expect_true(all(x$viz_data %in% viz))
})

test_that("DownSample.tviblindi resamples a sparse viz_data in lockstep", {
  data <- tiny_data(n = 40)
  x <- tviblindi(data = data, labels = list(pop = rep(c("a", "b"), 20)))
  viz <- Matrix::Matrix(tiny_data(n = 40, seed = 15), sparse = TRUE)
  Set_viz_data(x, viz)
  KNN(x, K = 10)
  DownSample(x, N = 15, K = 8)
  expect_equal(nrow(x$viz_data), nrow(x$data))
  expect_true(methods::is(x$viz_data, "Matrix"))
})
