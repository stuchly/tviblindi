#' Create downsampled version of tviblindi object
#'
#' This method takes a \code{tviblindi} object and creates a modified version with downsampled expression data.
#' 
#' However, working with the full data (without downsampling), if possible, is encouraged.
#'
#' @param tv \code{tviblindi}-class object
#' @param N integer: size of the domwnsampled data. Default value is \code{10000}
#' @param k integer: number of nearest neighbors to estimate density if density-based downsampling is used. Default value is \code{10}
#' @param method string: name of downsampling method to apply. Either \code{density} (for density-based downsampling) or \code{labels} (for downsampling biased against the removal of small labelled populations).
#' Default value is \code{density}
#' @param intrinsic_dim integer: expected intrinsic dimension of data. (Higher dimension under-estimates the density of sparse regions.) Default value is \code{2}
#' @param labels_name string: name of labels vector to use by method \code{labels}. Default value is '\code{default}'
#' @param threshold integer: threshold value for minimum population size in method \code{labels}. Default value is \code{floor(N / 100)}
#' @param ignore_labels string or string vector: names of labelled populations that should be left out entirely from sampling, for method \code{labels}. Default value is \code{NULL}
#'
#' @details
#' Density-based downsampling (for \code{method == 'density'}) should reduce over-abundant dense populations.
#' Label-based downsampling (for \code{method == 'labels'}) prevents the vanishing of smaller populations while downsampling, and works with a previously specified vector of population labels.
#'
#' @export
Downsample <- function(
  tv,
  N = 10000,
  k = 10,
  method = 'density',
  intrinsic_dim = 2,
  labels_name = 'default',
  threshold = floor(N / 100),
  ignore_labels = NULL
) {
  if (length(method) != 1 || !is.character(method))
    stop('Invalid format of "method"')
  if (!method[1] %in% c('density', 'labels'))
    stop(paste0('Invalid downsampling method: ', method))
  
  if (method == 'labels' && !labels_name %in% names(tv$labels))
    stop(paste0('Invalid labels vector namel: ', labels_name))
  
  max_N <- nrow(tv$data)
  if (N > max_N) {
    warning(paste0('Desired N is larger than number of vertices, set to ', max_N))
    N <- max_N
  }
  
  origin_labels <- if (!is.null(tv$origin)) tv$origin_label else NULL
  
  original_N <- nrow(tv$data)
  
  if (method[1] == 'density') {
    if (is.null(tv$kNN))
      stop('For "density" method of downsampling, k-NN matrix is needed')
    
    Vd <- pi^(intrinsic_dim / 2) / gamma(intrinsic_dim / 2 + 1)
    dens <- k / (max_N * Vd)
    dens <- dens / tv$kNN$dist[, k]^intrinsic_dim
    density_s <- sort(dens)
    cdf <- rev(cumsum(1.0 / rev(density_s)))
    boundary <- N / cdf[1]
    if (boundary > density_s[1]) {
      targets <- (N - 1:length(density_s)) / cdf
      boundary <- targets[which.min(targets - density_s > 0)]
    }
    s <- which(boundary / dens > stats::runif(length(dens)))
    .msg(paste0('Downsampled data to ', N, ' events using density'))
    
  } else if (method[1] == 'labels') {
    labels <- as.numeric(tv$labels[[labels_name]])
    if (!is.null(ignore_labels))
      labels[labels == labels[unique(as.character(tv$labels) %in% ignore_labels)]] <- -1
    labels.sorted <- as.numeric(names(sort(table(labels))))
    labels.sorted <- labels.sorted[labels.sorted != -1]
    
    if (N / length(labels.sorted) < threshold)
      .msg_alt_bad(paste0('N too small for the set threshold of ', threshold))
    
    labels.idcs <- lapply(labels.sorted, function(label) which(labels == label))
    
    out <- list()
    for (i in 1:length(labels.idcs)) {
      remaining_categories <- length(labels.idcs) - i
      remaining_N <- N - length(unlist(out))
      this_N <- length(labels.idcs[[i]])
      if (this_N < threshold) {
        out[[i]] <- labels.idcs[[i]]
      } else {
        selection_size <- min(remaining_N %/% (remaining_categories + 1), this_N)
        out[[i]] <- sample(x = labels.idcs[[i]], size = selection_size, replace = FALSE)
      }
    }
    
    .msg(paste0('Downsampled data to ', nrow(tv$data), ' events using label bias'))
    
    s <- unlist(out)
    if (length(s) < nrow(tv$data))
      .msg_alt_bad('Downsampled N lower than requested')
  }
  
  if (!is.null(origin_labels)) {
    .msg('Re-indexing cells-of-origin')
    
    origin_names <- names(tv$origin)
    
    tv$origin <- NULL
    for (idx.label in seq_len(origin_labels)) {
      label <- origin_labels[idx.label]
      name <- origin_names[idx.label]
      SetOrigin(tv, label = label, name = name)
    }
  }
  
  gc(verbose = FALSE)
  
  tv <- tviblindi(data = data, labels = labels)
  tv$origin <- origin
  tv$events_sel <- events_sel
  tv$downsampling <- downsampling
  tv$layout <- layout
  tv$labels <- labels
  tv
}
