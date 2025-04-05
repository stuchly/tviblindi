make_valid_fcs<-function(exprs,fcs=NULL,desc1=NULL){
  ##require(flowCore)
  if(!is.null(fcs)) fcs<-flowCore::flowFrame(exprs, parameters = flowCore::parameters(fcs)) else fcs<-flowCore::flowFrame(exprs)
##  for (i in 1:ncol(exprs)) if (min(exprs[,i])<0) exprs[,i]<-exprs[,i]-min(exprs[,i])
  params <- flowCore::parameters(fcs)
  pd <-NULL
  cols <- as.vector(pd$name)
  idxs <- match(cols, pd$name)

  if (any(is.na(idxs))) {
    stop("Invalid column specifier")
  }


  keyval <- list()
  for (channel_number in 1:ncol(exprs)){
    channel_name<-colnames(exprs)[channel_number]
    if (is.null(desc1)) desc1<-colnames(exprs)[channel_number]
    channel_id <- paste("$P", channel_number, sep = "")
    channel_range <- max(exprs[,channel_number]) + 1
    channel_min<-min(0,min(exprs[,channel_number])-1)
    plist <- matrix(c(channel_name, desc1[channel_number], channel_range,
                      channel_min, channel_range - 1))
    rownames(plist) <- c("name", "desc", "range", "minRange",
                         "maxRange")
    colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))




    keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
    keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
    keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
    keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
    keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
  }

  params@data <- data.frame(pd)
  if(!is.null(fcs)) fcs<-flowCore::flowFrame(exprs, parameters = params) else fcs<-flowCore::flowFrame(exprs,parameters = params)
  flowCore::keyword(fcs) <- keyval


  return(fcs)
}
