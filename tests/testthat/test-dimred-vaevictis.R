vaevictis_available <- function() {
  requireNamespace("reticulate", quietly = TRUE) &&
    isTRUE(tryCatch(reticulate::py_module_available("vaevictis"), error = function(e) FALSE))
}

tiny_tviblindi <- function() {
  set.seed(4)
  data <- matrix(rnorm(80 * 5), ncol = 5)
  labels <- sample(c("a", "b", "c"), 80, replace = TRUE)
  tviblindi(data = data, labels = list(pop = labels))
}

test_that("DimRed(method='vaevictis') trains and returns a finite [N,2] layout, and sets x$vae", {
  skip_if_not(vaevictis_available(), "vaevictis python module not available")

  x <- tiny_tviblindi()
  DimRed(
    x,
    method = "vaevictis",
    dim = 2,
    enc_shape = c(8, 8),
    dec_shape = c(8, 8),
    batch_size = 16L,
    epochs = 1L,
    patience = 1L,
    upsample = NULL,
    K = 5
  )

  got <- x$layout[[1]]
  expect_equal(dim(got), c(nrow(x$data), 2))
  expect_true(all(is.finite(got)))
  expect_false(is.null(x$vae))
})

test_that("load_model path does not set x$vae, matching pre-refactor behavior", {
  skip_if_not(vaevictis_available(), "vaevictis python module not available")
  skip("requires a pre-saved model fixture; documents the preserved quirk, not exercised here")
})
