SymmetriseMatrices.tviblindi <- function(
  tv,
  slots,
  name = NULL
) {
  for (slot in slots) {
    if (is.null(name)) {
      tv[[slot]] <- igraph::as_adjacency_matrix(
        igraph::graph_from_adjacency_matrix(
          tv[[slot]],
          weighted = TRUE, mode = 'max'
        ),
        attr = 'weight'
      )
    } else {
      tv[[slot]][[name]] <- igraph::as_adjacency_matrix(
        igraph::graph_from_adjacency_matrix(
          tv[[slot]][[name]],
          weighted = TRUE, mode = 'max'
        ),
        attr = 'weight'
      )
    }
  }
  invisible(tv)
}

SymmetriseMatrices <- function(tv, slots, name) UseMethod('SymmetriseMatrices', tv)

OrientMatrices.tviblindi <- function(
  tv,
  slots,
  name = NULL,
  pseudotime_name = NULL,
  base = NULL,
  breaks = NULL,
  suffix = ''
) {
  for (slot in slots) {
    
    ## Extract edges with weights
    N <- length(tv$pseudotime[[1]]$res)
    
    triplets <- Matrix::summary(
      if (is.null(name)) tv[[slot]] else tv[[slot]][[name]]
    )
    colnames(triplets) <- c('from', 'to', 'weight')
    
    ## Only keep those edges which lead to nodes with higher pseudotime value
    oriented <- triplets[which(tv$pseudotime[[pseudotime_name]]$res[triplets$from] < tv$pseudotime[[pseudotime_name]]$res[triplets$to]), ]
    
    if (!is.null(breaks)) {
      ## Apply penalty for short-circuiting
      b <- as.numeric(cut(as.numeric(as.factor(tv$pseudotime[[pseudotime_name]]$res)), breaks))
      b <- b[oriented$to] - b[oriented$from]
      oriented$weight <- oriented$weight / base^b
    }
    
    ## Graph Laplacian
    A  <- Matrix::sparseMatrix(
      i = oriented$from,
      j = oriented$to,
      x = oriented$weight,
      dims = c(N, N)
    )
    .d <- Matrix::Diagonal(x = (Matrix::rowSums(A))^(-1))
    oriented <- .d %*% A
    
    if (is.null(name)) {
      tv[[paste0(slot, suffix)]] <- oriented
    } else {
      tv[[slot]][[paste0(name, suffix)]] <- oriented
    }
  }
  invisible(tv)
}

OrientMatrices <- function(tv, slots, name, pseudotime_name, base, breaks, suffix) UseMethod('OrientMatrices', tv)

#' Compute transition model for expression data in a \code{tviblindi} object
#'
#' This method adds a new slot \code{trans} to an object of class \code{tviblindi}, containing transition probabilities associated with graph edges.
#' A chosen kernel is applied to compute transition probabilities from distances in \code{kNN}, with parameter \code{epsilon}.
#'
#' @param tv \code{tviblindi}-class object with \code{kNN}
#' @param name string: name of transition model to generate
#' @param origin_name string: name of cell-of-origin to use
#' @param k integer: number of nearest neighbours to use for calculating transition probabilities (upper bound is \code{k} of \code{tv$kNN})
#' @param kernel string: name of kernel function to use.
#' One of \code{Exp} (exponential), \code{SE} (standard error), \code{Lap} (Laplacian). Default value is \code{Exp}
#' @param epsilon numeric: kernel function parameter (epsilon). Default value is \code{NULL}, whereby a good epsilon is estimated automatically
#' 
#' @references
#' \insertRef{Bates2019}{tviblindi}
#'
#' @export
ModelTransitions.tviblindi <- function(
  tv,
  name,
  origin_name,
  k,
  kernel = 'Exp',
  epsilon = NULL
) {
  stopifnot(!is.null(tv$kNN))
  stopifnot(is.character(kernel) && kernel[1] %in% c('Exp', 'SE', 'Lap', 'ExpM', 'ExpMr', 'SEMr', 'ExpSigma'))
  
  N <- nrow(tv$kNN$idcs)
  
  .msg('Constructing distance matrix')
  tv$dist <- rbind(
    rep(1:N, each = k),
    as.numeric(t(tv$kNN$idcs[, 2:(k + 1)] + 1)),
    as.numeric(t(tv$kNN$dist[, 2:(k + 1)]))
  )
  dist <- tv$dist
  tv$dist <- Matrix::sparseMatrix(
    i = tv$dist[1, ],
    j = tv$dist[2, ],
    x = tv$dist[3, ],
    dims = c(N, N)
  )
  .msg('Symmetrising distance matrix')
  SymmetriseMatrices(tv, 'dist')
  gc(verbose = FALSE)
  
  if (is.null(tv$trans)) tv$trans <- list()
  
  .msg('Computing transition matrix')
  if (kernel[1] == 'SE') {
    if (is.null(epsilon)) epsilon <- 2 * median(dist[3, ])^2 # under-estimating
    tv$trans[[name]] <- dist
    tv$trans[[name]][3, ] <- exp(-tv$trans[[name]][3, ]^2 / epsilon)
  }
  if (kernel[1] == 'Lap') {
    if (is.null(epsilon)) epsilon <- median(dist[3, ]) # under-estimating
    tv$trans[[name]] <- dist
    tv$trans[[name]][3, ] <- exp(-tv$trans[3, ] / epsilon)
  }
  if (kernel[1] == 'Exp') {
    if (is.null(epsilon)) epsilon <- 2 * median(dist[3, ])^2 # under-estimating
    tv$trans[[name]] <- dist
    tv$trans[[name]][3, ] <- exp(-tv$trans[[name]][3, ] / epsilon)
  }
  if (kernel[1] == 'ExpM') {
    if (is.null(epsilon)) epsilon <- max(dist[3, ]) # under-estimating
    tv$trans[[name]] <- dist
    dist[3, ] <- exp(-dist[3, ] / epsilon)
  }
  if (kernel[1] == 'EmpMr') {
    d <- matrix(dist[3, ], nrow = N)
    m <- apply(d, MARGIN = 1, max)
    tv$trans[[name]] <- dist
    tv$trans[[name]][3, ] <- as.numeric(t(t(d) / m))
    tv$trans[[name]][3, ] <- exp(-tv$trans[[name]][3, ])
  }
  if (kernel[1] == 'SEMr') {
    d <- matrix(dist[3, ], nrow = N)
    m <- apply(d, MARGIN = 1, max)
    tv$trans[[name]] <- dist
    tv$trans[[name]][3, ] <- as.numeric(t(t(d^2) / m))
    tv$trans[[name]][3, ] <- exp(-tv$trans[[name]][3, ])
  }
  if (epsilon == 0) stop('Epsilon is zero')
  gc(verbose = FALSE)
  
  tv$trans[[name]] <- Matrix::sparseMatrix(
    i = tv$trans[[name]][1, ],
    j = tv$trans[[name]][2, ],
    x = tv$trans[[name]][3, ],
    dims = c(N, N)
  )
  gc(verbose = FALSE)
  
  .msg('Symmetrising transition matrix')
  SymmetriseMatrices(tv, 'trans', name)#, sym_mode = sym_mode)
  gc(verbose = FALSE)
  return(invisible(tv))
}

ModelTransitions <- function(tv, name, origin_name, k, kernel, epsilon) UseMethod('ModelTransitions', tv)

#' Assign distance-from-origin to each vertex for expression data in \code{tviblindi}
#'
#' This method adds a new slot \code{pseudotime} to an object of class \code{tviblindi}.
#'
#' Distance from cell-of-origin to each vertex is computed using an existing transition matrix.
#' Solving a system of linear equations, the average number of transitions needed to reach each vertex in \code{data} is computed.
#' If the number of vertices in \code{data} is less than a pre-set threshold, an analytical (exact) solution is computed; otherwise we apply conjugate gradient to solve numerically.
#'
#' @param tv \code{tviblindi}-class object with \code{kNN}, \code{dist} and \code{trans}
#' @param name string: name of transition model to use
#' @param origin_name string or integer: name of cell-of-origin to use
#' @param max_analytical integer: maximum number of vertices in data (without origin) for which we compute an exact solution. Default value is \code{1000}
#' @param nb_it integer: maximum number of iterations for conjugate gradient method. Default value is \code{1500}
#' @param tolerated_error number: error tolerance for conjugate gradient method. Default value is \code{1e-6}
#' 
#' @references
#' 
#' \insertRef{Bates2019}{tviblindi}
#' 
#' \insertRef{Bates2014}{tviblindi}
#' 
#' \insertRef{Eddelbuettel2013}{tviblindi}
#' 
#' @export
AssignDistances.tviblindi <- function(
  tv,
  name,
  origin_name,
  max_analytical = 1000,
  n_iter = 1500,
  tolerated_error = 1E-6
) {
  stopifnot(!is.null(tv$origin))
  
  ## Calculate pseudotime for each vertex
  
  N <- nrow(tv$trans[[name]])
  .D <- Matrix::rowSums(tv$trans[[name]])
  D <- Matrix::Diagonal(x = .D) # degree matrix D
  Dm <- Matrix::Diagonal(x = (.D)^(-1/2)) # normalised
  L  <- Dm %*% (D - tv$trans[[name]]) %*% Dm # symmetric normalised Laplacian
  rm(D)
  
  idcs.unknown <- which(!(1:nrow(L) %in% tv$origin[[origin_name]]))
  L <- L[idcs.unknown, idcs.unknown] # remove non-origin vertices
  
  B <- Matrix::rowSums((Dm %*% tv$trans[[name]]) * tv$dist) # D^{-1/2} %*% T 
  B <- matrix(B[idcs.unknown], nrow = length(idcs.unknown)) # column matrix without non-origin elements
  guess <- rep(0, N) # initial guess for solution
  
  if (length(idcs.unknown) > max_analytical) {
    .msg('Solving for pseudotime numerically')
    if (class(L) != 'dgCMatrix') stop('Use dgcMatrix')
    res <- solve_conjugate_gradient(
      L, B, guess[idcs.unknown], n_iter, tolerated_error
    )
  } else {
    message('Solving for pseudotime analytically')
    res <- list(
      x = Matrix::solve(L, B)[, 1],
      nb_it = 'exact',
      error = 'exact'
    )
  }
  pseudotime <- matrix(NA, nrow = N)
  pseudotime[tv$origin[[origin_name]], ] <- 0
  pseudotime[idcs.unknown, ] <- res$x
  pseudotime <- as.vector(Dm %*% pseudotime) # back to solution of standard Laplacian
  if (is.null(tv$pseudotime)) tv$pseudotime <- list()
  tv$pseudotime[[name]] <- list(
    res = pseudotime,
    nb_it = res$nb_it,
    error = res$error
  )
  invisible(tv)
}

AssignDistances <- function(tv, name, origin_name, max_analytical, n_iter, tolerated_error) UseMethod('AssignDistances', tv)

#' Infer pseudotime for expression data in \code{tviblindi}
#'
#' This method adds a new slot \code{pseudotime} to an object of class \code{tviblindi}.
#' A pseudotime value is assigned to each vertex in \code{data}, based on the cell-of-origin and transition probabilities between vertices.
#'
#' @param tv \code{tviblindi}-class object
#' @param name string: name of generated transition model. Default value is '\code{default}'
#' @param origin_name string: name of cell-of-origin to use. Default value is \code{1}
#' @param k integer: nearest neighbour count to use for computing transition probabilities. Upper bound is \code{k} of \code{tv$kNN}. Default value is \code{30}
#' @param kernel string: name of kernel function to use for constructing the transition model. One of \code{Exp} (exponential), \code{SE} (standard error), \code{Lap} (Laplacian). Default value is \code{Exp}
#' @param epsilon number: kernel function parameter (epsilon). Default value is \code{NULL}, whereby a good value of epsilon is estimated automatically
#' @param base numeric: penalty of jumping 'too far ahead' in pseudotime, to prevent 'short-circuiting'. Default value is \code{1.5}
#' @param breaks integer: number of bins for pseudotime values, for the subsequent simulation of random walks. Default value is \code{100}
#' @param n_iter maximum number of iterations for assignment of distance-from-origin to each vertex as the average number of transitions from origin. Usually solved using conjugate gradient (unless the dataset is small). Default value is \code{1500}
#'
#' @details
#' To avoid short-circuiting in subsequently simulated random walks (which model developmental trajectories), events are divided into equally sized bins (number specified by \code{breaks}).
#' The probability of jumping \code{k} bins ahead is then penalised by \code{base}^(-k).
#'
#' @references
#' \insertRef{Bates2019}{tviblindi}
#' 
#' \insertRef{Bates2014}{tviblindi}
#' 
#' \insertRef{Eddelbuettel2013}{tviblindi}
#' 
#' \insertRef{Csardi2006}{tviblindi}
#'
#' @export
ComputePseudotime.tviblindi <- function(
  tv,
  name = 'default',
  origin_name = 1,
  k = 30,
  kernel = 'Exp',
  epsilon = NULL,
  base = 1.5,
  breaks = 100,
  n_iter = 1500
) {
  
  if (length(tv$origin) == 0)
    stop('Cell-of-origin not set')
  if (is.null(tv$kNN))
    stop('kNN matrix not computed')
  
  max_k <- ncol(tv$kNN$idcs) - 1
  
  if (k > max_k) {
    k <- max_k
    message(paste0('k is larger than nearest neighbour count in kNN matrix, reducing k to ', max_k))
  }
  
  .msg('(1/3) Computing transition model')
  ModelTransitions(tv, name = name, origin_name = origin_name, k = k, kernel = kernel, epsilon = epsilon)#, sym_mode = sym_mode)
  gc(verbose = FALSE)
  .msg('(2/3) Computing pseudotime values')
  AssignDistances(tv, name = name, origin_name = origin_name, n_iter = n_iter)
  gc(verbose = FALSE)
  .msg('(3/3) Orienting transition model by pseudotime...')
  OrientMatrices(tv, 'dist', base = base, pseudotime_name = name, breaks = breaks)
  OrientMatrices(tv, 'trans', name = name, pseudotime_name = name, base = base, breaks = breaks)
  gc(verbose = FALSE)
  .msg_alt(paste0('Pseudotime computation error: ', tv$pseudotime[[name]]$error))
  
  if (is.null(tv$trans_kernel))
    tv$trans_kernel <- list()
  tv$trans_kernel[[name]] <- kernel
  
  if (is.null(tv$trans_name_origin))
    tv$trans_name_origin <- list()
  tv$trans_name_origin[[name]] <- origin_name
  
  invisible(tv)
}

ComputePseudotime <- function(tv, name, origin_name, k, kernel, epsilon, base, breaks, n_iter) UseMethod('ComputePseudotime', tv)