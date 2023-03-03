.onLoad<-function(libname,pckgname){
  if (reticulate::py_available("vaevictis")) reticulate::import("vaevictis")

  return(invisible(TRUE))
}
