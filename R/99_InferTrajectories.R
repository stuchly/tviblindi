#' Trajectory inference in single-cell expression data
#'
#' This function uses the \code{tviblindi} trajectory-inference pipeline to discover tentative developmental pathways in biological expression data (cytometry, scRNA-seq, CITE-seq ADT).
#' It returns an object of class \code{tviblindi}, containing the ordering of cells along a pseudotime and simulated random walks through the dataset.
#' To categorise walks by terminal nodes and classify them by differences in direction, use the function \code{Interactive} to launch an interactive session in your web browser.
#' This will allow you to look how different markers are expressed at each pseudotemporal segment and to export pathways of interest.
#'
#' ## Input data and annotations
#'
#' Input data are be provided in the form of a coordinate matrix.
#' You can either pass this matrix to \code{InferTrajectories} directly or, when working with cytometry data, you can give the path to an FCS file (and, optionally, row indices of events you want to use).
#' 
#' Other arguments to the \code{InferTrajectories} function include a vector of manually assigned population labels per cell and the label corresponding to the earliest progenitor cell population.
#' This is the population in which simulated developmental pathways always begin.
#' 
#' ## Model parameters
#' 
#' Furthermore, the list of parameters uncludes the neighbour count for a \code{k}-nearest-neighbours graph construction, number of denoising iterations to eliminate technical noise, size of SOM grid for clustering, parameters of a transition model, method to use for computing a 2-dimensional layout of the data and number of random walks to simulate.
#' For ease of use, these parameters have default values that give good results for some mass cytometry datasets.
#' 
#' To understand each of the constituent steps, set parameters more carefully and use multiple label vectors, transition models and 2-dimensional layouts, use the following functions to run each step separately:
#' 
#' * \code{tviblindi} for specifying input data and manually assigned labels,
#' 
#' * \code{SetOrigin} for specofying a population of origin,
#' 
#' * \code{Downsample} for downsampling your data if needed (better not to downsample whenever possible),
#' 
#' * \code{ConstructkNNG} for constructing a \code{k}-nearest-neighbours graps,
#' 
#' * \code{Denoise} for reducing technical variability in input data,
#' 
#' * \code{AddLayout} for reducing dimensionality of input data and producing a 2-dimensional visualisation,
#' 
#' * \code{Plot} for checking quality of the 2-dimensional layout,
#' 
#' * \code{Cluster} for creating a self-organising map (SOM) of data using \code{FlowSOM} prior to filtration,
#' 
#' * \code{Filter} for computing a witness-complex filtration of the high-dimensional data,
#' 
#' * \code{ComputePseudotime} for modelling transition probabilities between cells and computing pseudotime for each cell,
#' 
#' * \code{SimulateRandomWalks} for simulating stochastic walks governed by the transition model and pseudotemporal ordering,
#' 
#' ## Connectome analysis
#' 
#' In addition to a classical trajectory-inference analysis, you can let \code{tviblindi} cluster the expresion data using the Louvain method for community detection and then estimate flux in the network.
#' Use the function \code{Connectome} to do this.
#'
#' @param expression optional numeric matrix: coordinate matrix of expression data. Rows correspond to cells, columns correspond to markers. Columns must be named. Default value is \code{NULL}
#' @param fcs_path optional string: path to FCS file which contains the expression data. If \code{expression} is left as \code{NULL}, the expression matrix is extracted from this file. Default value is \code{NULL}
#' @param marker_names optional string vector: alternative expression matrix column names to the ones extracted from FCS file. Default value is \code{NULL}
#' @param fcs_rows optional integer vector: indices of expression matrix rows of FCS file that were taken (if \code{expression} is specified) or that should be taken (if \code{expression} is left as \code{NULL}). Default value is \code{NULL} (no subsetting)
#' @param labels_per_event string or factor vector: vector of population labels per each event (row of expression matrix)
#' @param origin_label string or integer: name population of origin or index of cell-of-origin
#' @param analysis_name string: name of \code{tviblindi} analysis. Default value is '*tviblindi TI analysis*'
#' @param n_neighbours_knn integer: \code{k} parameter for building \code{k}-NN graph. Default value is \code{100}
#' @param n_neighbours_denoise integer: number of nearest neighbours for each point to use for data denoising (if denoising is used). Maximum \code{n_neighbours_knn}. Default value is \code{30}
#' @param n_iter_denoise integer: number of denoising iterations. Set to \code{0} to disable denoising. Default value is \code{1}
#' @param somgrid_width integer: width of SOM grid for \code{FlowSOM} clustering of expression data. Default value is \code{25}
#' @param somgrid_width integer: height of SOM grid for \code{FlowSOM} clustering of expression data. Default value is \code{25}
#' @param transition_model_kernel string: name of kernel function to use for computing transition probabilities between cells.
#' One of \code{Exp} (exponential), \code{SE} (standard error), \code{Lap} (Laplacian). Default value is \code{Exp}
#' @param n_neighbours_transition_model integer: integer: number of nearest neighbours for each point to use for transition model building. Maximum \code{n_neighbours_knn}. Default value is \code{30}
#' @param layout_method string: name of dimension-reduction method for producing a 2-dimensional layout. One of \code{vaevictis}, \code{umap} and \code{tsne}. \code{vaevictis} is recommended, but it requires the \code{vaevictis} Python package and \code{reticulate} to be installed prior to calling \code{InferTrajectories}. Default value is \code{umap}
#' @param n_walks integer: number of random walks to simulate. Default value is \code{5000}
#'
#' @export
InferTrajectories <- function(
  expression                    = NULL,
  fcs_path                      = NULL,
  marker_names                  = NULL,
  fcs_rows                      = NULL,
  labels_per_event,
  origin_label,
  analysis_name                 = 'tviblindi TI analysis',
  n_neighbours_knn              = 100,
  n_neighbours_denoise          = 30,
  n_iter_denoise                = 1,
  somgrid_width                 = 25,
  somgrid_height                = 25,
  transition_model_kernel       = 'Exp',
  n_neighbours_transition_model = 30,
  layout_method                 = 'umap',
  n_walks                       = 5000,
  verbose                       = TRUE
) {
  
  if (!is.null(expression) && (!is.numeric(expression) || !is.matrix(expression) || is.null(colnames(expression)) || nrow(expression) < 2 || ncol(expression) < 2))
    stop('"expression" must be a numeric matrix with multiple rows, multiple columns and column name set')
  if (!is.null(fcs_path) && is.character(fcs_path) && is.atomic(fcs_path) && (!file.exists(fcs_path) || toUpper(substr(fcs_path, nchar(fcs_path) - 3, nchar(fcs_path))) != '.FCS'))
    stop('"fcs_path" must be a valid path to an FCS file')
  
  if (!is.null(fcs_path) && is.null(expression)) {
    expression <- flowCore::read.FCS(fcs_path)@exprs
    if (!is.null(fcs_rows)) {
      if (!all(fcs_rows) %in% seq_len(nrow(expression)) && any(duplicated(fcs_rows)))
        stop('"fcs_rows" should be a vector of unique indices in the range from 1 to number of rows of the "expression" matrix')
      expression <- expression[fcs_rows, ]
    }
  }
  
  if (!is.null(marker_names) && (!is.character(marker_names) || length(marker_names) != ncol(expression)))
    stop('"marker_names" must be a string vector of length equal to number of cÏolumns of the "expression" matrix')
  
  if (!is.null(marker_names))  
    colnames(expression) <- marker_names
  
  if (length(labels_per_event) != nrow(expression))
    stop('Length of labels vector must be equal to the number of rows of the expression data')
  if (is.numeric(origin_label) && (origin_label > nrow(expression)))
    stop('If "origin_label" is numeric, it must be a valid row index of the expression data')
  if (is.character(origin_label) && !origin_label %in% unique(labels_per_event))
    stop('If "origin_label" is a string, it must be one of the labels from "labels_per_event"')
  
  if (verbose)
    .msg_name(analysis_name)
  
  tv <- tviblindi(
    data            = expression,
    labels          = labels_per_event,
    fcs_path        = fcs_path,
    fcs_subset_idcs = fcs_rows,
    analysis_name   = analysis_name
  )
  
  if (verbose)
    .msg('Setting cell-of-origin')
  suppressMessages(SetOrigin(tv, label = origin_label))
  
  if (verbose)
    .msg('Constructing k-NN graph')
  suppressMessages(ConstructkNNG(tv, k = n_neighbours_knn))
  
  if (n_iter_denoise > 0) {
    if (verbose)
      .msg('Denoising expression data')
    suppressMessages(Denoise(tv, k = n_neighbours_denoise, n_iter = n_iter_denoise))
  }
  
  if (verbose)
    .msg('Creating 2-d layout')
  suppressMessages(AddLayout(tv, method = layout_method))
  
  if (verbose)
    .msg('Clustering')
  suppressMessages(Cluster(tv, xdim = somgrid_width, ydim = somgrid_height))
  
  if (verbose)
    .msg('Computing filtration')
  suppressMessages(Filter(tv))
  
  if (verbose)
    .msg('Computing pseudotemporal ordering of cells')
  suppressMessages(ComputePseudotime(tv, k = n_neighbours_transition_model, kernel = transition_model_kernel))
  
  if (verbose)
    .msg('Simulating random walks')
  suppressMessages(SimulateRandomWalks(tv, n_walks = n_walks))
  
  if (verbose) {
    .msg_alt_good('Done')
    .msg_alt('Returning a tviblindi-analysis object')
    .msg_alt('(To explore developmental pathways in your expression data, use the function "Interactive")')
  }
  
  tv
}