#' Dimensionality reduction via vaevictis (deep autoencoder)
#'
#' Built-in "vaevictis" \code{DimRed} method. Registered automatically at package load; see
#' \code{\link{register_dimred_method}} for the method-plugin contract this follows.
#'
#' @param x tviblindi class object; may be mutated (\code{x$vae} is set to the trained/loaded model).
#' @param data numeric matrix; the resolved input data (raw or denoised), as passed by \code{DimRed}.
#' @param dim integer (default 2); dimension of reduced data.
#' @param vsplit double (default 0.1); percentage of data used as validation step.
#' @param enc_shape integer vector (default \code{c(128,128,128)}); shape (depth and width) of the encoder.
#' @param dec_shape integer vector (default \code{c(128,128,128)}); shape (depth and width) of decoder.
#' @param perplexity double (default 10.); perplexity for tsne regularisation, see https://www.nature.com/articles/s41467-018-04368-5.
#' @param batch_size integer (default 512); training batch size.
#' @param epochs integer (default 100); maximum number of training epochs.
#' @param patience integer (default 1); early-stopping patience.
#' @param ivis_pretrain integer (default 0); ivis pretraining epochs.
#' @param ww numeric vector (default \code{c(10,10,1,1)}); loss weights - tsne regularisation, ivis pn loss, reconstruction error, KL divergence.
#' @param margin double (default 1); ivis pn loss margin.
#' @param shuffle logical (default FALSE); shuffle data before validation split; involves recomputation of the KNN matrix.
#' @param K integer (default 30); number of nearest neighbors for the KNN matrix used in training.
#' @param load_model character vector of 2 components (optional); paths to files created by \code{x$vae$save(file1,file2)} - model is loaded and applied instead of training a new one.
#' @param upsample named list \code{list(N=5000,cluster=15,method="kmeans")} or \code{NULL}; sample events by clusters (involves recomputation of the KNN matrix); if \code{NULL} nothing happens. \code{N} events per cluster, \code{cluster} number of clusters, \code{method} clustering method (either "CLARA" or "kmeans").
#' @param labels_name name of the label vector to use if \code{upsample$cluster} is \code{NULL} and \code{x} has multiple label vectors.
#' @param ... unused; absorbs any extra arguments DimRed's dispatch may forward.
#'
#' @return numeric matrix of reduced coordinates.
#'
#' @export
dimred_vaevictis <- function(x, data, dim = 2, vsplit = 0.1,
                              enc_shape = c(128, 128, 128),
                              dec_shape = c(128, 128, 128),
                              perplexity = 10., batch_size = 512L,
                              epochs = 100L, patience = 1L,
                              ivis_pretrain = 0,
                              ww = c(10., 10., 1., 1.), margin = 1.,
                              shuffle = FALSE, K = 30, load_model = NULL,
                              upsample = list(N = 5000, cluster = 15, method = "kmeans"),
                              labels_name = names(x$labels)[1], ...) {
  vv <- get_vaevictis()

  if (!is.null(load_model)) {
    model <- vv$loadModel(config_file = load_model[1], weights_file = load_model[2])
    layout <- model[[2]](data)
    # NOTE: x$vae is intentionally left untouched here, matching the pre-refactor
    # DimRed.tviblindi behavior (the load_model path never persisted the loaded
    # model onto x$vae, unlike the trained-model path below).
    return(layout)
  }

  if (!is.null(upsample)) {
    if (is.null(upsample$method)) upsample$method <- "CLARA"
    if (is.null(upsample$cluster)) {
      labl <- x$labels[[labels_name]]
    } else {
      if (upsample$method == "kmeans") {
        message("running kmeans clustering...")
        labl <- kmeans(data, centers = upsample$cluster)$cluster
        message("~done\n")
      } else {
        message("running clara clustering...")
        labl <- cluster::clara(data, k = upsample$cluster)$clustering
        message("~done\n")
      }
    }
    ss <- .upsample.labels(labl, N = upsample$N, takeall = upsample$takeall)
    knn_loc <- KNN.annoy(data[ss, ], K, 150)$IND
    layout <- vv$dimred(
      data[ss, ], as.integer(dim), vsplit, as.integer(enc_shape),
      as.integer(dec_shape), perplexity, as.integer(batch_size),
      as.integer(epochs), as.integer(patience), as.integer(ivis_pretrain),
      ww, "euclidean", margin, K, knn_loc
    )
  } else {
    if (shuffle) sshuf <- sample(nrow(data)) else sshuf <- 1:nrow(data)
    if (shuffle) {
      knn.plc <- KNN.annoy(data[sshuf, ], K, 150)$IND
    } else {
      if (!is.null(x$KNN)) knn.plc <- KofRawN(x$KNN, K) else knn.plc <- KNN.annoy(data[sshuf, ], K, 150)$IND
    }
    layout <- vv$dimred(
      data[sshuf, ], as.integer(dim), vsplit, as.integer(enc_shape),
      as.integer(dec_shape), perplexity, as.integer(batch_size),
      as.integer(epochs), as.integer(patience), as.integer(ivis_pretrain),
      ww, "euclidean", margin, K, knn.plc
    )
  }
  x$vae <- layout[[3]]
  layout[[2]](data)
}
