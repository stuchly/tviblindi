tiny_tviblindi <- function() {
  set.seed(1)
  data <- matrix(rnorm(40 * 3), ncol = 3)
  labels <- sample(c("a", "b"), 40, replace = TRUE)
  tviblindi(data = data, labels = list(pop = labels))
}

test_that("built-in methods are registered at load time", {
  expect_true(all(c("vaevictis", "diffuse", "umap") %in% list_dimred_methods()))
})

test_that("register_dimred_method + DimRed dispatch round-trips through the registry", {
  register_dimred_method(".test_identity", function(x, data, ...) {
    x$.test_marker <- "was here"
    cbind(data[, 1], data[, 2])
  }, overwrite = TRUE)

  x <- tiny_tviblindi()
  DimRed(x, method = ".test_identity")

  expect_equal(names(x$layout), "1_.test_identity")
  expect_equal(dim(x$layout[[1]]), c(nrow(x$data), 2))
  expect_equal(x$.test_marker, "was here") # method mutated x directly, as documented
})

test_that("registering an existing name without overwrite errors", {
  register_dimred_method(".test_dup", function(x, data, ...) data, overwrite = TRUE)
  expect_error(
    register_dimred_method(".test_dup", function(x, data, ...) data, overwrite = FALSE),
    "already registered"
  )
  # overwrite=TRUE succeeds
  expect_silent(register_dimred_method(".test_dup", function(x, data, ...) data, overwrite = TRUE))
})

test_that("unknown method gives an informative error listing available methods", {
  x <- tiny_tviblindi()
  expect_error(DimRed(x, method = "definitely_not_a_method"), "Unknown dimred method.*vaevictis")
})

test_that("manual layout= bypasses method dispatch entirely", {
  register_dimred_method(".test_should_not_run", function(x, data, ...) stop("must not be called"), overwrite = TRUE)
  x <- tiny_tviblindi()
  manual <- matrix(1:(nrow(x$data) * 2), ncol = 2)
  DimRed(x, method = ".test_should_not_run", layout = manual)
  expect_equal(x$layout[[1]], manual)
})

test_that("successive DimRed calls append and index layouts correctly", {
  register_dimred_method(".test_identity2", function(x, data, ...) data[, 1:2], overwrite = TRUE)
  x <- tiny_tviblindi()
  DimRed(x, method = ".test_identity2")
  DimRed(x, method = ".test_identity2")
  expect_equal(names(x$layout), c("1_.test_identity2", "2_.test_identity2"))
})

test_that("method-specific arguments are forwarded via ...", {
  register_dimred_method(".test_args", function(x, data, scale = 1, ...) data[, 1:2] * scale, overwrite = TRUE)
  x <- tiny_tviblindi()
  DimRed(x, method = ".test_args", scale = 10)
  expect_equal(x$layout[[1]], x$data[, 1:2] * 10)
})
