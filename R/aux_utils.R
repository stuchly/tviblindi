.upsample.labels<-function(labels,N=10000,replace=TRUE,takeall="ungated"){
    out<-NULL
    out<-which(labels %in% takeall)
    for (i in unique(labels[!(labels %in% takeall)])) out<-c(out,sample(which(labels==i),size=ifelse(replace,N,min(N,length(which(labels==i)))),replace=replace))
    return(out[sample(length(out))])
}
