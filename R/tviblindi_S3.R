tviblindi<-function(data,labels){
    new_tviblindi(data,labels)
}

new_tviblindi<-function(data,labels,keep.intermediate=FALSE){
    stopifnot(is.matrix(data))
    stopifnot(is.vector(labels) && length(labels)==nrow(data) && (is.factor(labels) || is.character(labels)))
    stopifnot(is.logical(keep.intermediate))

    out<-new.env(hash=TRUE)
    out$data=data
    out$labels<-labels
    out$keep<-keep.intermediate
    out$pseudotime<-NULL
    out$filtration<-NULL
    out$boundary<-NULL
    out$walks<-NULL
    out$KNN<-NULL

    structure(out,class="tviblindi")
}
