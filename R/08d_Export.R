#' Add additional channels to FCS file with info about pathways of interest
#'
#' \code{add_pathways_info_to_FCS} requires a \code{flowFrame} object, pseudotime values, all simulated random walks, indices of trajectories of interest,
#' indices of selected events (in case the analysed data correspond only to a subset of data in the \code{flowFrame}) and a a numeric or character ID for the 
#' specified collection of trajectories.
#' 
#' A valid \code{flowFrame} object with additional channels is returned.
#' Firstly, the channel \code{'ID_which_event'} contains a binary values tagging which events are included in the trajectories of interest (1000) and which were not (-1).
#' Secondly, the channel \code{'ID_local_pseudotime'} contains pseudotime values of vertices comprising the trajectories of interest is created.
#' In both cases \code{'ID'} is replaced with the corresponding parameter value.
#' You do not need to call this function directly.
#'
#' @param ff \code{flowFrame} object (loaded using package \code{flowCore}) corresponding to FCS with the analysed data
#' @param pseudotime list: must contain a numeric vector in the slot \code{res}, corresponding to a pseudotime value for each event included in the analysis
#' @param all_walks list of random walks, consisting of node indices in slot \code{v} and start indices referring to entries from \code{v} in \code{starts}
#' @param trajectories_of_interest integer vector: walk indices, with entries from \code{1:length(walks$starts)}
#' @param fcs_subset_idcs integer vector (optional): indices of selected events
#' Those are the row numbers of expression matrix from the FCS file which correspond to events included in \code{data}
#' Default value is \code{1:nrow(data)}, which means all events were selected
#' @param id numeric or character: name of the batch of trajectory vertices being appended
#'
#' @return \code{flowFrame} object (corresponds to an enhanced FCS file)
#'
#' @references
#' \insertRef{Ellis2019}{tviblindi}
#'
#' @export
add_pathways_info_to_FCS <- function(
  ff,
  pseudotime,
  all_walks,
  trajectories_of_interest,
  fcs_subset_idcs = NULL,
  id
) {
  pp            <- unique(select_walks_nodes(all_walks, trajectories_of_interest))
  out           <- matrix(-1, nrow = nrow(ff@exprs), ncol = 2)
  colnames(out) <- c(paste(id, 'which_event',      sep = '_'),
                     paste(id, 'local_pseudotime', sep = '_'))
  if (is.null(fcs_subset_idcs)) fcs_subset_idcs <- 1:nrow(ff@exprs)
  
  out[fcs_subset_idcs[pp], 1] <- 1000
  out[fcs_subset_idcs[pp], 2] <- as.numeric(as.factor(pseudotime$res[pp]))
  
  make_valid_fcs(
    cbind(ff@exprs, out),
    desc1 = as.character(ff@parameters@data$desc)
  )
}


fcs.add_col <- function(
  ff,
  new_col,
  colname = 'label'
) {
  if (class(ff) != 'flowFrame') stop('Parameter "ff" is not a valid flowFrame object')
  
  efcs <- ff@exprs
  N    <- nrow(efcs)
  len  <- length(new_col)
  
  if (N != len)  stop(paste0('Number of rows of expression matrix is ', N, ', whereas length of new column is ', len, '.'))
  
  params <- ff@parameters
  pd     <- flowCore::pData(params)
  cols   <- as.vector(pd$name)
  idcs   <- match(cols, pd$name)
  
  if (any(is.na(idcs))) stop('Invalid column specifier')
  
  channel_number     <- ncol(ff) + 1
  channel_id         <- paste0('$P', channel_number)
  channel_name       <- colname
  channel_range      <- max(new_col) + 1
  channel_min        <- min(0, min(new_col) - 1)
  plist              <- matrix(c(channel_name, channel_name, channel_range,
                                 channel_min, channel_range - 1))
  rownames(plist)    <- c('name', 'desc', 'range', 'minRange', 'maxRange')
  colnames(plist)    <- c(channel_id)
  pd                 <- rbind(pd, t(plist))
  flowCore::pData(params) <- pd
  channel_names      <- colnames(efcs)
  efcs.mod           <- cbind(efcs, new_col)
  colnames(efcs.mod) <- c(channel_names, colname)
  
  ff.mod             <- flowCore::flowFrame(efcs.mod, params, description = flowCore::description(ff))
  
  keyval                                      <- list()
  keyval[[paste0('$P', channel_number, 'B')]] <- '32'
  keyval[[paste0('$P', channel_number, 'R')]] <- toString(channel_range)
  keyval[[paste0('$P', channel_number, 'E')]] <- '0,0'
  keyval[[paste0('$P', channel_number, 'N')]] <- channel_name
  keyval[[paste0('$P', channel_number, 'S')]] <- channel_name
  flowCore::keyword(ff.mod)                             <- keyval
  
  flowCoreP_Rmax <- paste0('flowCore_$P', channel_number, 'Rmax')
  flowCoreP_Rmin <- paste0('flowCore_$P', channel_number, 'Rmin')
  
  flowCore::description(ff.mod)[flowCoreP_Rmax] <- max(20000, flowCore::description(ff.mod)$`flowCore_$P1Rmax`)
  flowCore::description(ff.mod)[flowCoreP_Rmin] <- 0
  
  ff.mod
}

make_valid_fcs <- function(
  exprs,
  fcs = NULL,
  desc1 = NULL
) {
  fcs    <- if (!is.null(fcs)) { flowCore::flowFrame(exprs, parameters = flowCore::parameters(fcs)) } else { flowCore::flowFrame(exprs) }
  params <- flowCore::parameters(fcs)
  pd     <- NULL
  cols   <- as.vector(pd$name)
  idxs   <- match(cols, pd$name)
  
  if (any(is.na(idxs))) {
    stop('Invalid column specifier')
  }
  keyval <- list()
  for (channel_number in 1:ncol(exprs)){
    channel_name<-colnames(exprs)[channel_number]  
    if (is.null(desc1)) desc1<-colnames(exprs)[channel_number] 
    channel_id <- paste('$P', channel_number, sep = '')
    channel_range <- max(exprs[,channel_number]) + 1
    channel_min<-min(0,min(exprs[,channel_number])-1)
    plist <- matrix(c(channel_name, desc1[channel_number], channel_range, 
                      channel_min, channel_range - 1))
    rownames(plist) <- c('name', 'desc', 'range', 'minRange', 
                         'maxRange')
    colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))

    keyval[[paste('$P', channel_number, 'B', sep = '')]] <- '32'
    keyval[[paste('$P', channel_number, 'R', sep = '')]] <- toString(channel_range)
    keyval[[paste('$P', channel_number, 'E', sep = '')]] <- '0,0'
    keyval[[paste('$P', channel_number, 'N', sep = '')]] <- channel_name
    keyval[[paste('$P', channel_number, 'S', sep = '')]] <- channel_name
  }
  
  params@data           <- data.frame(pd)
  if(!is.null(fcs)) fcs <- flowCore::flowFrame(exprs, parameters = params) else fcs <- flowCore::flowFrame(exprs,parameters = params)
  flowCore::keyword(fcs)          <- keyval

  fcs
}
