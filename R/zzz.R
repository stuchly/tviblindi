.onLoad<-function(libname,pckgname){
  if (any(c(!reticulate::py_module_available("tensorflow"),
            !reticulate::py_module_available("numba"),
            !reticulate::py_module_available("numpy"),
            !reticulate::py_module_available("vaevictis")))){
    message("Not all required python modules are installed. Install following modules?")
    if (!reticulate::py_module_available("tensorflow")) message("tensorflow")
    if (!reticulate::py_module_available("numba")) message("numba")
    if (!reticulate::py_module_available("numpy")) message("numpy")
    if (!reticulate::py_module_available("vaevictis")) message("vaevictis - dimensional reduction for tviblindi")
    yn<-readline("Proceed? (Yes/no)")
    if (yn!="Yes"){
      warning("function DimRed functionality missing!")
      return(invisible(FALSE))
    }
    if (!reticulate::py_module_available("tensorflow")) reticulate::py_install("tensorflow",pip=TRUE)
    if (!reticulate::py_module_available("numba")) reticulate::py_install("numba",pip=TRUE)
    if (!reticulate::py_module_available("numpy")) reticulate::py_install("numpy",pip=TRUE)
    if (!reticulate::py_module_available("vaevictis")) reticulate::py_install("git+https://github.com/stuchly/vaevictis.git@master",pip=TRUE)
  }
  
  return(invisible(TRUE))
}