launch_shiny_old <- function(input_fcs_path,
                         origin,
                         pseudotime,
                         coords,
                         coords_clusters,
                         clusters,
                         filtration,
                         layout,
                         walks_raw,
                         b,
                         rb,
                         event_sel = NULL,
                         export = TRUE,
                         overwrite = FALSE) {
  if (export) {
    shiny_inputs_dir <- "tviblindi_tmp"
    if (dir.exists(shiny_inputs_dir) & !overwrite) {
      message("Shiny inputs folder exists. Proceed?")
      readline()
    }
    
    dir.create(shiny_inputs_dir)
    message("Exporting objects for Shiny")
    message("-> path to analysed FCS file")
    saveRDS(input_fcs_path, file.path(shiny_inputs_dir, "input_fcs_path.RDS"))
    message("-> cell of origin index")
    saveRDS(origin, file.path(shiny_inputs_dir, "origin.RDS"))
    message("-> pseudotime values")
    saveRDS(pseudotime, file.path(shiny_inputs_dir, "pseudotime.RDS"))
    message("-> expression matrix")
    saveRDS(coords, file.path(shiny_inputs_dir, "coords.RDS"))
    message("-> clustered expression matrix")
    saveRDS(coords_clusters, file.path(shiny_inputs_dir, "coords_clusters.RDS"))
    message("-> cluster definitions")
    saveRDS(clusters, file.path(shiny_inputs_dir, "clusters.RDS"))
    message("-> filtration")
    saveRDS(filtration, file.path(shiny_inputs_dir, "filtration.RDS"))
    message("-> 2-dimensional layout")
    saveRDS(layout, file.path(shiny_inputs_dir, "layout.RDS"))
    message("-> simulated random walks")
    saveRDS(walks_raw, file.path(shiny_inputs_dir, "walks_raw.RDS"))
    message("-> boundary matrix")
    saveRDS(b, file.path(shiny_inputs_dir, "b.RDS"))
    message("-> reduced boundary matrix")
    saveRDS(rb, file.path(shiny_inputs_dir, "rb.RDS"))
    message("-> indices of events selected for analysis")
    saveRDS(event_sel, file.path(shiny_inputs_dir, "event_sel.RDS"))
  }
  
  app <- shiny::shinyApp(shiny_ui, shiny_server)
  shiny::runApp(app)
}

launch_shiny <- function(input_fcs_path,
                         tviblindi_s3,
                         event_sel = NULL,
                         export = TRUE,
                         overwrite = TRUE) {
  
  tmp_folder <- 'tviblindi_tmp'
  tv_path    <- file.path(tmp_folder, 'tv.RDS')
  tviblindi_s3_name <- deparse(substitute(tviblindi_s3))
  message(paste0("WARNING: your tviblindi analysis object '", tviblindi_s3_name, "' will be deleted from this environment and saved as an .RDS file '", tv_path, "'.\nAfter the Shiny session ends, load your analysis object using function 'readRDS'"))
  
  tviblindi_s3$event_sel <- event_sel
  
  if (export) {
    shiny_inputs_dir <- "tviblindi_tmp"
    if (dir.exists(shiny_inputs_dir) & !overwrite) {
      message("Shiny inputs folder exists. Proceed?")
      readline()
    }
  }
  
  suppressWarnings(dir.create(tmp_folder))
  
  message("Created 'tviblindi_tmp' data directory for Shiny UI")
  message(paste0('--> exporting analysis data to ', file.path(tmp_folder, 'tv.RDS')))
  saveRDS(input_fcs_path, file.path(tmp_folder, 'input_fcs_path.RDS'))
  saveRDS(tviblindi_s3,   tv_path)
  
  message('--> removing analysis data from parent environment')
  rm(list = tviblindi_s3_name, envir = sys.frame(-1))
  
  message('--> launching Shiny UI')
  app <- shiny::shinyApp(shiny_ui, shiny_server)
  shiny::runApp(app)
  
  
}

