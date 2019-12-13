addPathInfot2Fcs<-function(fcs,pseudotime,walks,walks.selection,event_sel=NULL,ID=NULL){
    pp<-unique(select_paths_points(walks,walks.selection))
    out<-matrix(-1,nrow=nrow(exprs(fcs)),ncol=2)
    colnames(out)<-c(paste(ID,"which_event",sep="_"),paste(ID,"local_pseudotime",sep="_"))
    if (is.null(event_sel)) event_sel<-1:nrow(exprs(fcs))

    out[event_sel[pp],1]<-1000
    out[event_sel[pp],2]<-as.numeric(as.factor(pseudotime$res[pp]))
    
    make_valid_fcs(cbind(exprs(fcs),out),desc1=as.character(fcs@parameters@data$desc))
}
