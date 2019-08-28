Xvalues <- function(colX, sorted, n_windows){
  a <- rep( NA , n_windows)
  for (i in 1:n_windows){
    r <- which(sorted==i)

    a[i] <- median(colX[r])
  }

  return(a)
}

plot_data_sets<-function(data,n_windows,pseudotime,paths,cls,clssel=unique(cls)){

}

sort_dots5<-function(pptag01, n_windows){
  probs<-(seq(0,1, 1/n_windows))^2
  probs<-probs/max(probs)

  br <- unique(as.numeric(quantile(pptag01,probs)))

  sorted <- as.numeric(cut(pptag01, breaks= br,include.lowest = TRUE))
  return(sorted)

}

plot_data<-function (data, n_windows, pp, parameters, tit = NULL, plines = NULL,
                     gate = NULL, add = FALSE, lty = 1)
{
    ## palette(c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5"))
    palette(brewer.pal(n=length(parameters),name="Set2"))
  pptag01 <- pp/max(pp)
  sorted <- sort_dots5(pptag01, n_windows)

  if (!is.null(gate)){
    mm<-data.frame(gate=gate,window=sorted)
    tt<-as.matrix(table(mm))
    ts<-colSums(tt)
    tt<-t(t(tt)/ts)
    N<-ncol(tt)
    xx<-barplot(tt,col=1:length(levels((gate))),ylim=c(-0.4,1))
    tt<-as.numeric(table(sorted))
    text(x=xx,y=-0.05,as.character(tt),cex=0.5,col="red",srt=90)
    legend("bottomleft", legend = levels((gate)),
           pch = 15,
           col = 1:length(levels((gate))),
           cex = 0.6)
  }
  CC <- 1:length(parameters)
  x1 <- 1:n_windows
  x <- x1/n_windows - 1/n_windows
  y <- Xvalues(data[, parameters[1]], sorted, n_windows)
  M <- max(y, na.rm = TRUE)
  m <- min(y, na.rm = TRUE)
  for (i in 1:length(parameters)) {
    k <- parameters[i]
    y <- Xvalues(data[, k], sorted, n_windows)
    if (M < max(y, na.rm = TRUE)) {
      M <- max(y, na.rm = TRUE)
    }
    if (m > min(y, na.rm = TRUE)) {
      m <- min(y, na.rm = TRUE)
    }
  }
  if (!add) {
    par(mar = c(4.5, 1, 4.1, 1))
    plot(x, y, type = "l", col = CC[1], ylim = c(m - 0.2 *
                                                   (M - m), 1.2 * M), main = tit, xlab = "", ylab = "",
         yaxt = "n", xlim = c(0, 1.1))
    if (FALSE) {
      pch <- as.numeric((gate))

      vv1 <- 0.65/n_windows
      for (ii in unique(sorted)) {
        ss <- which(sorted == ii)
        xv <- runif(length(ss)) * vv1
        yv <- runif(length(ss)) * 1/5 * M
        points(x[ii] + xv, y[ii] + yv, pch = pch[ss],
               cex = 0.5, col = pch[ss])
      }
      ## print(levels(gate))
      ## print(head(gate))
      ## print(as.numeric(head(gate)))
      legend("bottomleft", legend = levels((gate)),
             pch = 1:length(levels((gate))),
             col = 1:length(levels((gate))),
             cex = 0.5)
    }
  }
  legendP <- NULL
  for (i in 1:length(parameters)) {
    k <- parameters[i]
    y <- Xvalues(data[, k], sorted, n_windows)
    lines(x, y, col = CC[i], lty = lty)
    legendP <- c(legendP, colnames(data)[k])
  }
  colp <- c("green", "green", "blue", "blue", "red", "red")
  if (!is.null(plines))
    for (i in 1:length(plines)) abline(v = pptag01[plines[i]],
                                       lty = 2, col = colp[i], lwd = 0.5)
  legend("topleft", legend = legendP, pch = 15, col = CC, cex = 0.5)
  return(sorted)
}


## plot_data <- function(data, n_windows, pp, parameters,tit=NULL,plines=NULL,gate=NULL,add=FALSE,lty=1){   ## parameters = vector containing the numbers of the columns to plot
##   pptag01 <- pp /max(pp)
##   ##w <- create_windows(n_windows)
##   sorted <- sort_dots3(pptag01, n_windows)                     ## Sorting the points into the windows

##   ## CC <- rainbow(length(parameters))
##   CC<-1:length(parameters)
##   x1 <- 1:n_windows
##   x <- x1/n_windows - 1/n_windows
##   y <- Xvalues(data[,parameters[1]], sorted, n_windows)        ## Getting the values from the initial table

##   ## Computing the limits of the y axis
##   M <- max(y,na.rm=TRUE)
##   ## YY<<-y
##   ## print(M)
##   m <- min(y,na.rm=TRUE)
##   for(i in 1:length(parameters)){
##     k<-parameters[i]
##     y <- Xvalues(data[,k], sorted, n_windows)
##     ## y<-y/max(y)
##     ## print(max(y))
##     if(M < max(y,na.rm=TRUE)) {
##       M <- max(y,na.rm=TRUE)
##     }
##     if(m > min(y,na.rm=TRUE)){
##       m <- min(y,na.rm=TRUE)
##     }

##   }
##   ## plot the data
##   if (!add){
##       par(mar=c(4.5, 1, 4.1, 1))
##       plot(x,y,type="l",col=CC[1], ylim=c(m-0.2*(M-m),1.2*M),main=tit,xlab = "",ylab="",yaxt = "n",xlim=c(0,1.1))
##       if (!is.null(gate)){
##           pch<-as.numeric(as.factor(gate))
##           vv1<-0.65/n_windows
##           for (ii in unique(sorted)){
##               ss<-which(sorted==ii)
##               xv<-runif(length(ss))*vv1
##               yv<-runif(length(ss))*1/5*M
##               points(x[ii]+xv,y[ii]+yv,pch=pch[ss],cex=0.5,col=pch[ss])
##           }
##           legend("bottomleft",legend=levels(as.factor(gate)),pch=as.numeric(as.factor(levels(as.factor(gate)))),col=as.numeric(as.factor(levels(as.factor(gate)))),cex=0.5)
##       }
##   }
##   legendP <- NULL
##   for (i in 1:length(parameters)){
##     k<-parameters[i]
##     y <- Xvalues(data[,k], sorted, n_windows)
##     ## y<-y/max(y)
##     lines(x, y, col = CC[i],lty=lty)
##     legendP <- c(legendP, colnames(data)[k])
##     ##print(k)
##   }
##   colp<-c("green","green","blue","blue","red","red")
##   if (!is.null(plines)) for (i in 1:length(plines)) abline(v=pptag01[plines[i]],lty=2,col=colp[i],lwd=0.5)
##   legend("topleft", legend = legendP, pch = 15, col=CC,cex=0.5)
##   return(sorted)
## }





create_windows <- function(n_windows,quanti=FALSE){
    if (!quanti){
        windows <- list()
        k1 <- 0:n_windows
        k1 <- k/n_windows
        for (i in 1:n_windows) {

            windows[[i]] <- c(k1[i],k1[i+1])
        }
    } else {

    }
    return(windows)
}


sort_dots3 <- function(pptag01, n_windows){
    br <- seq(0,1, 1/n_windows)

    sorted <- as.numeric(cut(pptag01, breaks= br,include.lowest = TRUE))
    return(sorted)

}

sort_dots4 <- function(pptag01, n_windows){
    br <- as.numeric(quantile(pptag01,seq(0,1, 1/n_windows)))
    sorted <- as.numeric(cut(pptag01, breaks= br))
    return(sorted)

}
