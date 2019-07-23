## TDA_scan<-function(data,time,breaks=20,maxdimension=1){
##     print(max(time))
##     if (length(breaks)==1) breaks<- as.numeric(quantile(time,seq(0,1,length.out = breaks)))
##     ##if (length(breaks)==1) breaks<- seq(0,max(time),length.out=breaks)
##     out<-list()[1:(length(breaks)-1)]
##     for (i in length(breaks):2){
##         print(i)
##         print(breaks[i])
##         ## data_loc<-data[time>=breaks[i-1] & time<breaks[i],]
##          data_loc<-data[ time>=breaks[i],]
##         ##maxscale<-max(knn.adj.raw(data_loc,1)$DIST)

##         ## FltRips <- ripsFiltration(X = data_loc, maxdimension = maxdimension,
##         ##                   maxscale = maxscale, dist = "euclidean", library = "GUDHI",
##         ##                   printProgress = TRUE)
##         FiltAlp<-alphaComplexFiltration(X=data_loc,printProgress = TRUE)
##         out[[i-1]]<- filtrationDiag(filtration = FiltAlp, maxdimension = maxdimension,library = "Dionysus", location = TRUE, printProgress = TRUE)
##     }
##     return(out)
## }
TDA_scan<-function(data,time,breaks=20,maxdimension=1){
  print(max(time))
  if (length(breaks)==1) breaks<- as.numeric(quantile(time,seq(0,1,length.out = breaks)))
  ##if (length(breaks)==1) breaks<- seq(0,max(time),length.out=breaks)
  data_loc_out<-out<-list()[1:(length(breaks)-2)]
  for (i in (length(breaks)-1):2){
    print(i)
    ## data_loc<-data[time>=breaks[i-1] & time<breaks[i],]
    data_loc<-data[ time>=breaks[i],]
    data_loc_out[[i-1]]<-data_loc

    ##maxscale<-max(knn.adj.raw(data_loc,1)$DIST)

    ## FltRips <- ripsFiltration(X = data_loc, maxdimension = maxdimension,
    ##                   maxscale = maxscale, dist = "euclidean", library = "GUDHI",
    ##                   printProgress = TRUE)
    FiltAlp<-alphaComplexFiltration(X=data_loc,printProgress = TRUE)
    out[[i-1]]<- filtrationDiag(filtration = FiltAlp, maxdimension = maxdimension,library = "Dionysus", location = TRUE, printProgress = TRUE)
  }
  return(list(out,data_loc_out))
}
