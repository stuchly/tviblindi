.dimred_registry <- new.env(parent = emptyenv())

#' Register a dimensionality-reduction method for use with DimRed()
#'
#' \code{register_dimred_method} adds (or replaces) an entry in tviblindi's dimensionality-reduction
#' method registry, making it available as \code{DimRed(x, method=name, ...)}. Built-in methods
#' ("vaevictis", "diffuse", "umap") are registered through this exact same mechanism at package
#' load time, so third-party methods are first-class, not a special case.
#'
#' @param name character; the method name used as \code{DimRed}'s \code{method} argument.
#' @param fn function with signature \code{function(x, data, ...)}. \code{x} is the tviblindi
#' object (an environment - the function may set additional fields on it, e.g. to persist a
#' trained model, the way the built-in "vaevictis" method sets \code{x$vae}); \code{data} is the
#' resolved input matrix (raw or denoised, per \code{DimRed}'s \code{use.denoised}); \code{...}
#' are the method's own parameters, declared with their own defaults on \code{fn}'s formals.
#' Must return a numeric matrix of the reduced coordinates.
#' @param overwrite logical (default FALSE); if FALSE, registering an already-used name errors
#' instead of silently replacing it.
#'
#' @return invisible \code{NULL}.
#'
#' @export
register_dimred_method <- function(name, fn, overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1)
  stopifnot(is.function(fn))
  if (!overwrite && exists(name, envir = .dimred_registry, inherits = FALSE)) {
    stop(sprintf(
      "dimred method '%s' is already registered; pass overwrite=TRUE to replace it",
      name
    ))
  }
  assign(name, fn, envir = .dimred_registry)
  invisible(NULL)
}

#' List registered dimensionality-reduction methods
#'
#' @return character vector of method names usable as \code{DimRed(x, method=...)}.
#'
#' @export
list_dimred_methods <- function() {
  sort(ls(envir = .dimred_registry))
}

.get_dimred_method <- function(name) {
  if (!exists(name, envir = .dimred_registry, inherits = FALSE)) {
    stop(sprintf(
      "Unknown dimred method '%s'. Available methods: %s",
      name, paste(list_dimred_methods(), collapse = ", ")
    ))
  }
  get(name, envir = .dimred_registry)
}
