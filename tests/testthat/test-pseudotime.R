tiny_data <- function() {
  set.seed(5)
  matrix(rnorm(80 * 4), ncol = 4)
}

test_that("pseudotime() on a data matrix matches Pseudotime.tviblindi() on the same data", {
  data <- tiny_data()
  root <- 1L

  x <- tviblindi(data = data, labels = list(pop = rep("a", nrow(data))))
  Set_origin(x, label = root, origin_name = "root")
  KNN(x, K = 10)
  Pseudotime(x, K = 10, origin_name = "root", weighted = TRUE)
  expected <- x$pseudotime[["root"]]$res

  got <- pseudotime(data, root = root, K = 10, weighted = TRUE)

  expect_equal(as.numeric(got), as.numeric(expected))
})

test_that("root and terminal_cells can be given as rownames, resolved via rownames(data)", {
  data <- tiny_data()
  rownames(data) <- paste0("cell", seq_len(nrow(data)))

  by_name <- pseudotime(data, root = "cell1", terminal_cells = "cell5", K = 10)
  by_index <- pseudotime(data, root = 1L, terminal_cells = 5L, K = 10)

  expect_equal(as.numeric(by_name), as.numeric(by_index))
  expect_equal(names(by_name), rownames(data))
})

test_that("unknown cell name errors clearly instead of silently misindexing", {
  data <- tiny_data()
  rownames(data) <- paste0("cell", seq_len(nrow(data)))
  expect_error(pseudotime(data, root = "not_a_cell", K = 10), "not_a_cell")
})

test_that("multiple terminal_cells match a direct ground-truth harmonic solve", {
  # regression test for the graph_assign.R boundary-value fixes made earlier this session:
  # with >=2 targets, both the correctness of unlabeled-cell pseudotime and the targets'
  # own reported values depend on those fixes (a namespace-qualification bug in the fix's
  # Matrix::t() usage was also caught and fixed by this test, see graph_assign.R). Verified
  # against the raw (unnormalized) hitting-distance recursion directly, not a heuristic
  # inequality - on a tiny graph a huge terminal value legitimately dominates the whole
  # field, so "interior values stay small" is not a safe assumption to test against.
  data <- tiny_data()
  root <- 1L
  terminal_cells <- c(5L, 10L, 15L)
  terminal_values <- c(20, 40, 60)

  knn <- KNN.annoy(data, 10, 150)
  d <- knn.raw2adj(knn)
  dsym <- knn.spadj2sym(knn.adj2spadj(d))
  sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = "SEMer"))
  Deg <- Matrix::rowSums(sim)
  P <- Matrix::Diagonal(x = 1 / Deg) %*% sim

  unlabeled <- setdiff(seq_len(nrow(data)), c(root, terminal_cells))
  I <- Matrix::Diagonal(nrow(data))
  LHS <- (I - P)[unlabeled, unlabeled]
  # weighted hitting distance: local cost c_i = sum_j P_ij * w_ij, using dsym as edge weights
  c_i <- Matrix::rowSums((Matrix::Diagonal(x = 1 / Deg) %*% sim) * dsym)[unlabeled]
  RHS <- c_i + as.numeric(P[unlabeled, terminal_cells, drop = FALSE] %*% terminal_values)
  ground_truth <- as.numeric(Matrix::solve(as.matrix(LHS), RHS))

  got <- pseudotime(data, root = root, terminal_cells = terminal_cells,
                     terminal_values = terminal_values, K = 10)

  expect_equal(unname(got[terminal_cells]), terminal_values)
  expect_equal(got[[root]], 0)
  expect_equal(unname(got[unlabeled]), ground_truth, tolerance = 1e-6)
})

test_that("a precomputed graph= is used as-is instead of rebuilding one from data", {
  data <- tiny_data()
  rownames(data) <- paste0("cell", seq_len(nrow(data)))

  knn <- KNN.annoy(data, 10, 150)
  d <- knn.raw2adj(knn)
  sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = "SEMer"))

  from_graph <- pseudotime(data = data, root = 1L, graph = sim, weighted = FALSE)
  from_data <- pseudotime(data = data, root = 1L, K = 10, weighted = FALSE)

  expect_equal(as.numeric(from_graph), as.numeric(from_data))
})

test_that("solver diagnostics are attached without polluting the primary numeric return", {
  data <- tiny_data()
  got <- pseudotime(data, root = 1L, K = 10)

  expect_true(is.numeric(got))
  expect_false(is.list(got))
  diag <- attr(got, "diagnostics")
  expect_true(all(c("nb_it", "error") %in% names(diag)))
})
