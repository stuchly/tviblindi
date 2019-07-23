normalize.perc<-function(data,id=rep(1,nrow(data)),ctrl=unique(id),probs=0.99){
    ss<-which(id %in% ctrl)
    qp<-apply(data[ss,],quantile,MARGIN=2,probs=probs)
    data<-t(t(data)/qp)
}
