#' Run  \emph{catGenes} under a Shiny Application
#'
#' @author Domingos Cardoso
#'
#' @description This function runs the \emph{catGenes} application. It will open
#' the application by using the default internet browser, or the viewer if running
#' from RStudio.
#'
#' @note If you have questions about the application, please email directly to
#' Domingos Cardoso at \email{cardosobot@@gmail.com}.
#'
#' @param launch.browser Run catGenes in browser
#' @param port If launch.browser is FALSE, specify port to run catGenes
#'
#' @examples if (interactive()) {
#' run_catGenes()
#' }
#'
#'
#' @export

run_catGenes <- function(launch.browser = TRUE, port = getOption("shiny.port")) {
  appDir <- system.file("catGeneshiny", package = "catGenes")
  if (appDir == "") {
    stop("Could not find catGeneshiny directory. Try re-installing `catGenes`.", call. = FALSE)
  }

  return(shiny::runApp(appDir, launch.browser = launch.browser, port = port))
}
