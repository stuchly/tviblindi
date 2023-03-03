.onLoad<-function(libname,pckgname){
  vae<-reticulate::import("vaevictis")

  return(invisible(TRUE))
}
