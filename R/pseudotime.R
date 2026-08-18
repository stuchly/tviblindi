#' Compute pseudotime (expected hitting time/distance) on a KNN graph
#'
#' \code{pseudotime} computes the expected hitting time (\code{weighted=FALSE}) or expected
#' hitting distance (\code{weighted=TRUE}) from a root cell to every other cell, over a
#' K-nearest-neighbor graph - the same graph-Laplacian boundary-value solver used internally
#' by \code{\link{Pseudotime.tviblindi}}, exposed here as a standalone function that only needs
#' a plain data matrix (or a precomputed graph), not a \code{tviblindi} object. Optionally, one
#' or more terminal cells can be pinned to a fixed (typically large) pseudotime value - mirroring
#' how e.g. Palantir's \code{early_cell}/\code{terminal_states} work - to mark known endpoints.
#'
#' @param data numeric matrix (cells x features); required unless \code{graph} is supplied.
#' Used both to build the KNN graph (unless \code{graph} is given) and to resolve character
#' \code{root}/\code{terminal_cells} against \code{rownames(data)}.
#' @param root integer index or character rowname of the root (origin) cell.
#' @param terminal_cells integer indices or character rownames of terminal cells (optional).
#' Pseudotime at these cells is pinned to \code{terminal_values} rather than solved for.
#' @param terminal_values numeric; pseudotime value(s) to pin \code{terminal_cells} to. Defaults
#' to \code{1e9} for every terminal cell if not supplied.
#' @param weighted logical (default TRUE); expected hitting distance (edge-weighted) if TRUE,
#' expected hitting time (step count) if FALSE.
#' @param K integer (default 30); number of nearest neighbors.
#' @param kernel character (default "SEMer"); similarity kernel, see \code{\link{knn.adj2spadjsim}}.
#' @param symmetrize character (default "max"); how to symmetrize the KNN graph - one of
#' "max", "mean", "prob", "min", "none". See \code{\link{knn.spadj2sym}} and related.
#' @param graph similarity \code{sparseMatrix} (optional); a precomputed graph to use instead of
#' building one from \code{data} - e.g. to reuse the same graph across multiple roots. When
#' supplied, \code{K}/\code{kernel}/\code{symmetrize} are ignored.
#' @param weights matrix of edge weights (optional), used only when \code{graph} is supplied and
#' \code{weighted=TRUE}; see \code{\link{assign_distance}}.
#' @param method character (default "cg"); numerical solver, "cg" or "minres".
#' @param nb_it integer (default 1500); maximum solver iterations.
#' @param eps double (default 1e-15); solver error tolerance.
#' @param iguess numeric vector (optional); initial guess for the iterative solver.
#'
#' @return a numeric vector of pseudotime values, one per row of \code{data} (or per vertex of
#' \code{graph}), named by \code{rownames(data)} if present. Solver diagnostics (\code{nb_it},
#' \code{error}) are attached as \code{attr(result, "diagnostics")}.
#'
#' @export
pseudotime <- function(data = NULL,
                        root,
                        terminal_cells = NULL,
                        terminal_values = NULL,
                        weighted = TRUE,
                        K = 30,
                        kernel = "SEMer",
                        symmetrize = "max",
                        graph = NULL,
                        weights = NULL,
                        method = "cg",
                        nb_it = 1500,
                        eps = 1e-15,
                        iguess = NULL) {
  stopifnot(!is.null(data) || !is.null(graph))

  resolve <- function(cell) {
    if (is.character(cell)) {
      idx <- match(cell, rownames(data))
      if (anyNA(idx)) stop("cell name(s) not found in rownames(data): ", paste(cell[is.na(idx)], collapse = ", "))
      idx
    } else {
      as.integer(cell)
    }
  }
  root <- resolve(root)
  if (!is.null(terminal_cells)) terminal_cells <- resolve(terminal_cells)
  if (!is.null(terminal_cells) && is.null(terminal_values)) {
    terminal_values <- rep(1e9, length(terminal_cells))
  }

  if (is.null(graph)) {
    knn <- KNN.annoy(data, K, 150)
    d <- knn.raw2adj(knn)
    dsym <- knn.spadj2sym(knn.adj2spadj(d))
    sym_input <- TRUE
    if (symmetrize == "none") {
      sim <- knn.adj2spadjsim(d, kernel = kernel)
      sym_input <- FALSE
    } else if (symmetrize == "mean") {
      sim <- knn.spadj.symmetrize(knn.adj2spadjsim(d, kernel = kernel))
    } else if (symmetrize == "prob") {
      sim <- knn.spadj.symmetrize.P(knn.adj2spadjsim(d, kernel = kernel))
    } else if (symmetrize == "max") {
      sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = kernel))
    } else if (symmetrize == "min") {
      d2 <- t(summary(dsym))
      sim <- knn.adj2spadjsim1(d2, kernel = kernel)
      sym_input <- FALSE
    } else {
      stop("symmetrize must be one of 'max', 'mean', 'prob', 'min', 'none'")
    }
    weights <- if (weighted) dsym else NULL
  } else {
    sim <- graph
    sym_input <- Matrix::isSymmetric(sim)
    if (!weighted) weights <- NULL
  }

  res <- assign_distance(
    sim, root,
    weights = weights, nb_it = nb_it, iguess = iguess, eps = eps,
    sym = sym_input, method = method,
    target = terminal_cells, target_values = terminal_values
  )

  out <- as.numeric(res$res)
  if (!is.null(data) && !is.null(rownames(data))) names(out) <- rownames(data)
  attr(out, "diagnostics") <- list(nb_it = res$nb_it, error = res$error)
  out
}
