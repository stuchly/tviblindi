.onLoad <- function(libname, pkgname) {
  register_dimred_method("vaevictis", dimred_vaevictis, overwrite = TRUE)
  register_dimred_method("diffuse", dimred_diffuse, overwrite = TRUE)
  register_dimred_method("umap", dimred_umap, overwrite = TRUE)
  register_dimred_method("phate", dimred_phate, overwrite = TRUE)
}
