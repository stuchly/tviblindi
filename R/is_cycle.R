is_cycle<-function(x,bb,rbb){
    x<-unlist(bb$cmplx[unlist(rbb$boundary[unlist(x)])])
    return(x)
    return(all(as.vector(table(x))%%2==0))
}








