#' Add a column to a flowFrame object
#'
#' \code{add1column} appends a column of data to the \code{exprs} data matrix of a \code{flowFrame} object and returns the extended \code{flowFrame} object.
#' @param fcs \code{flowFrame} object.
#' @param col numeric vector of values to be appended.
#' @param colname name of the appended column.
#' @param exprs_in \code{exprs} data matrix to be extended. If NULL, \code{exprs} is extracted from \code{fcs}.
#' @return An extended \code{flowFrame} object.
#' @export

add1column <- function(fcs, col, colname="new data", exprs_in=NULL)
{
  if (is.null(exprs_in)) in_data <- exprs(fcs) else in_data <- exprs_in
  if (nrow(in_data)!=length(col)) stop("args do not match in size")
  params <- parameters(fcs)
  pd <- pData(params)
  cols <- as.vector(pd$name)
  idxs <- match(cols, pd$name)

  if (any(is.na(idxs))) {
    stop("Invalid column specifier")
  }

  channel_number <- ncol(fcs) + 1
  channel_id <- paste("$P", channel_number, sep = "")
  channel_name <- colname
  channel_range <- max(col) + 1
  plist <- matrix(c(channel_name, channel_name, channel_range,
                    0, channel_range - 1))
  rownames(plist) <- c("name", "desc", "range", "minRange",
                       "maxRange")
  colnames(plist) <- c(channel_id)
  pd <- rbind(pd, t(plist))
  pData(params) <- pd
  cnames<-colnames(in_data)
  out_data <- cbind(in_data, col)
  colnames(out_data)<-c(cnames,colname)
  out_frame <- flowFrame(out_data, params, description = description(fcs))
  keyval <- list()
  keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
  keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
  keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
  keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
  keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
  keyword(out_frame) <- keyval
  return(out_frame)
}
