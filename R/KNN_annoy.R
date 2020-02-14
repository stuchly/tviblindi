KNN.annoy<-function(data,K,trees=150){
  res<-knn_annoy(data,K,trees)
  res<-list(IND=do.call("rbind",res[[1]]),DIST=do.call("rbind",res[[2]]))
  return(res)
}