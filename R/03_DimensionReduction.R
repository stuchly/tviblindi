#' Apply dimension reduction to produce a 2-d layout of expression data
#'
#' @param tv \code{tviblindi}-class object
#' @param method string: dimension-reduction method. One of \code{vaevictis} (encoder with t-SNE regularisation), \code{umap} or \code{tsne}. Default value is \code{vaevictis}
#' @param name string: name of layout to be generated. Default value is \code{default}
#' @param layout optional numeric matrix with number of rows equal to \code{nrow(tv$data)}: pre-computed layout matrix, to be added instead of creating a new layout. Default value is \code{NULL}
#' @param vsplit number: portion of data used as validation set in \code{vaevictis} method. Default value is \code{0.1}
#' @param encoder_shape integer vector of length 2: shape (depth and width) of encoder in \code{vaevictis}. Default value is \code{c(128,128,128)}
#' @param decoder_shape integer vector of length 2: shape (depth and width) of decoder in \code{vaevictis}. Default value is \code{c(128,128,128)}
#' @param perplexity number: perplexity parameter for t-SNE regularisation in \code{vaevictis} or in t-SNE. Default value is \code{10}
#' @param batch_size integer: batch size for \code{vaevictis} training. Default value is \code{512}
#' @param n_epochs integer: maximum number of epochs for training a \code{vaevictis} or \code{umap} model. Default is \code{100} for \code{vaevictis} and 200 for \code{umap}
#' @param n_neighbours integer: neighbour count for \code{umap}. Default value is \code{50}
#' @param patience integer: maximum patience parameter for \code{vaevictis} training (early stopping). Default value is \code{0}
#' @param alpha number: penalty parameter for t-SNE regularisation in \code{vaevictis}. Default value is \code{10}
#' @param load_model character vector of length 2: paths to files created by \code{tv$vae$save(file1, file2)} for loading a \code{vaevictis} model instead of training a new one.
#' Default value is \code{NULL}, which causes a new model to be trained
#'
#' @references
#' \insertRef{Nl2008}{tviblindi}
#' 
#' \insertRef{McInnes2018}{tviblindi}
#' 
#' \insertRef{Ushey2020}{tviblindi}
#' 
#' @export
AddLayout.tviblindi <- function(
  tv,
  method = 'vaevictis',
  name = 'default',
  layout = NULL,
  n_epochs = if (method == 'vaevictis') 100L else if (method == 'umap') 200 else NULL,
  n_neighbours = if (is.null(tv$kNN)) 50L else ncol(tv$kNN$idcs) - 1,
  load_model = NULL
) {
  if (!is.character(method) || (!method[1] %in% c('vaevictis', 'umap', 'tsne')))
    stop(paste0('Invalid dimension-reduction method: ', method))
  
  if (!is.null(layout)) {
    if (is.null(tv$layout)) tv$layout <- list()
    tv$layout[[name]] <- layout
    .msg('Appended existing layout to tviblindi object')
    return(invisible(tv))
  }
  
  if (method[1] == 'vaevictis') {
    if (is.null(tv$kNN))
      stop('k-NN matrix needs to be computed first for vaevictis')
    available <- function(pkg) grepl('TRUE', invisible(utils::capture.output(reticulate::py_module_available(pkg))))
    
    av.tf <- available('tensorflow')
    av.numba <- available('numba')
    av.numpy <- available('numpy')
    av.vae <- available('vaevictis')
    
    if ((sum(!c(av.tf, av.numba, av.numpy, av.vae)) -> n_missing) > 0) {
      .msg_alt_bad(paste0('The following ', n_missing, ' Python module', if (n_missing > 1) { 's' } else { '' }, ' required by vaevictis ',
        if (n_missing > 1) { 'are' } else { 'is' }, ' missing:\n',
          if (!av.tf)    { '     tensorflow\n' },
          if (!av.numba) { '     numba\n' },
          if (!av.numpy) { '     numpy\n' },
          if (!av.vae)   { '     vaevictis (git+https://github.com/davnovak/vaevictis.git@master)\n'}))
      yn <- readline("To install requirements via reticulate, type 'Yes'...")
      if (yn != 'Yes') {
        .msg_alt_bad('vaevictis dimension-reduction functionality missing')
        return(invisible(FALSE))
      }
      
      if (!av.tf)
        reticulate::py_install('tensorflow', pip = TRUE)
      if (!av.numba)
        reticulate::py_install('numba', pip = TRUE)
      if (!av.numpy)
        reticulate::py_install('numpy', pip = TRUE)
      if (!av.vae)
        reticulate::py_install('git+https://github.com/davnovak/vaevictis.git@master', pip = TRUE)
      
      .msg_alt_good('All required modules installed')
    }
    
    vv <- reticulate::import('vaevictis')
    if (!is.null(load_model)) {
      model <- vv$loadModel(config_file = load_model[1], weights_file = load_model[2])
      tv$layout <- model$encoder$callp(tv$data)$numpy()
      tv$vae <- model
    } else {
      layout <- vv$dimred(
        tv$data,
        epochs = as.integer(n_epochs),
        knn = tv$kNN$IND
      )
      if (is.null(tv$vae))    tv$vae <- list()
      if (is.null(tv$layout)) tv$layout <- list()
      tv$vae[[name]]    <- layout[[3]]
      tv$layout[[name]] <- layout[[1]]
    }
  } else if (method[1] == 'umap') {
    if (is.null(tv$layout)) tv$layout <- list()
    tv$layout[[name]] <- uwot::umap(tv$data, n_neighbors = n_neighbours, n_components = 2L)
  } else if (method[1] == 'tsne') {
    if (is.null(tv$layout)) tv$layout <- list()
    tv$layout[[name]] <- Rtsne::tsne(tv$data, k = 2)
  }
  tv$layout_method <- method
  gc(verbose = FALSE)
  
  invisible(tv)
}

AddLayout <- function(tv, method, name, layout, n_epochs, n_neighbours, load_model) UseMethod('AddLayout', tv)