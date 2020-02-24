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

        tv                  <- list()
        tv$origin           <- origin
        tv$pseudotime       <- pseudotime
        tv$data             <- coords
        tv$codes            <- coords_clusters
        tv$clusters         <- clusters
        tv$filtration       <- filtration
        tv$layout           <- layout
        tv$walks            <- walks_raw
        tv$boundary         <- b
        tv$reduced_boundary <- rb
        tv$event_sel        <- event_sel


        dir.create(shiny_inputs_dir)
        message("Exporting analysis object for Shiny")
        saveRDS(tv, file.path(shiny_inputs_dir, 'tv.RDS'))
        saveRDS(input_fcs_path, file.path(shiny_inputs_dir, 'input_fcs_path.RDS'))

    }

    app <- shiny::shinyApp(shiny_ui, shiny_server)
    shiny::runApp(app)
}

#' Launches shiny interface, modifies x$ShowAllFates
#'
#' \code{Filtration}
#' @param tviblindi_s3 tviblindi class object.
#' @param ShowAllFates bool; plot all theoretical fates, see \code{walks}
#' 
#'
#' @export
launch_shiny <- function(tviblindi_s3,ShowAllFates=TRUE) {
    
    tviblindi_s3$ShowAllFates<-ShowAllFates
    if (class(tviblindi_s3) != 'tviblindi') stop('Invalid tviblindi S3 object')
    tmp_folder        <- 'tviblindi_tmp'
    tviblindi_s3_name <- deparse(substitute(tviblindi_s3))
    suppressWarnings(dir.create(tmp_folder))
    tv_path           <- file.path(tmp_folder, 'tv.RDS')
    saveRDS(tviblindi_s3_name, tv_path)

    app <- shiny::shinyApp(shiny_ui, shiny_server)
    shiny::runApp(app)
}

