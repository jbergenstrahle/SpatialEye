#' @export
run <- function(dir) {
  shiny::addResourcePath("tiles", dir)
  shiny::runApp(appDir = system.file("", package = "spatialeye"))
}
