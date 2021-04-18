#' Launch interactive Shiny interface for tviblindi
#'
#' WARNING: running this function creates a folder names \code{tviblindi_tmp} in the working directory.
#'
#' @param tviblindi_s3 \code{tviblindi} class object
#' @param show_all_fates logical: whether all potential terminal nodes (according to the transition model) should be shown, regardless of whether any simulated walks actually terminate in them
#'
#' @references
#' 
#' \insertRef{Chang2019}{tviblindi}
#' 
#' \insertRef{Pederson2019}{tviblindi}
#' 
#' \insertRef{Perrier2020}{tviblindi}
#' 
#' \insertRef{Bailey2015}{tviblindi}
#' 
#' \insertRef{DeVries2016}{tviblindi}
#' 
#' \insertRef{Galili2015}{tviblindi}
#' 
#' \insertRef{Wickham2019}{tviblindi}
#' 
#' \insertRef{Ellis2019}{tviblindi}
#' 
#' \insertRef{Wickham2019a}{tviblindi}
#' 
#' \insertRef{Warnes2020}{tviblindi}
#' 
#' \insertRef{Kratochvil2020}{tviblindi}
#'
#' @export
Interactive <- function(
  tviblindi_s3,
  show_all_fates = TRUE
) {
  if (class(tviblindi_s3) != 'tviblindi')
    stop('Invalid tviblindi S3 object')
  if (is.null(tviblindi_s3$layout))
    stop('2-dimensional layout of data needed must be computed before launching the Shiny app')
  if (is.null(tviblindi_s3$filtration))
    stop('Filtration needs to be computed before launching the Shiny app')
  if (is.null(tviblindi_s3$pseudotime))
    stop('Pseudotime values need to be computed before launching the Shiny app')
  if (is.null(tviblindi_s3$walks))
    stop('Random walks need to be simulated before launching the Shiny app')

  tv$show_all_fates <- show_all_fates
  
  tmp_folder <- 'tviblindi_tmp'
  tviblindi_s3_name <- deparse(substitute(tviblindi_s3))
  suppressWarnings(dir.create(tmp_folder))
  tv_path <- file.path(tmp_folder, 'tv.RDS')
  saveRDS(tviblindi_s3_name, tv_path)
  
  app <- shiny::shinyApp(shiny_ui, shiny_server)
  shiny::runApp(app)
}

