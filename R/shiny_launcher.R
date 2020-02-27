#' Launches shiny interface
#'
#' \code{Filtration}
#' @param tviblindi_s3 tviblindi class object.
#' 
#' @details See \code{ToggleShowFates}.
#'
#' @export
launch_shiny <- function(tviblindi_s3, ShowAllFates = TRUE) {
    
    if (class(tviblindi_s3) != 'tviblindi') stop('Invalid tviblindi S3 object')
    tmp_folder        <- 'tviblindi_tmp'
    tviblindi_s3_name <- deparse(substitute(tviblindi_s3))
    suppressWarnings(dir.create(tmp_folder))
    tv_path           <- file.path(tmp_folder, 'tv.RDS')
    saveRDS(tviblindi_s3_name, tv_path)
    
    app <- shiny::shinyApp(shiny_ui, shiny_server)
    shiny::runApp(app)
}
