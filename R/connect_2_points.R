connect2points<-function(X,i0,i1,by=0.1){
    ## X<-cbind(X,0)
    L<-dim(X)[2]
    t<-seq(by,1,by=by)
    perp<-t*10
    handle_01<-cbind(rbind(matrix(rep(X[i0,],length(perp)),length(perp),byrow = TRUE),matrix(rep(X[i1,],length(perp)),length(perp),byrow = TRUE)),rep(perp,2))


    ##matrix(rep(X[i0,],length(t)),length(t),byrow = TRUE)+matrix(rep(X[i1,]-X[i0,],length(t)),length(t),byrow=TRUE)*t
    handle_01<-rbind(handle_01,cbind(matrix(rep(X[i0,],length(t)),length(t),byrow = TRUE)+matrix(rep(X[i1,]-X[i0,],length(t)),length(t),byrow=TRUE)*t,max(perp)))

   return(rbind(cbind(X,0),handle_01))
}
