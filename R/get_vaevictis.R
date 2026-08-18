get_vaevictis <- local({
  vaevictis_module <- NULL
  function() {
    if (is.null(vaevictis_module)) {
      if (!reticulate::py_available(initialize = FALSE)) {
        ## stop("Python not initialized. Please call reticulate::import('vaevictis') manually before calling this.")
      }
      vaevictis_module <- reticulate::import("vaevictis")
    }
    vaevictis_module
  }
})
