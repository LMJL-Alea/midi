#' Runs the MIDI Shiny web application
#'
#' This is a helper function to run the MIDI Shiny web application in the
#' default web browser.
#'
#' @return Nothing, but launches the Shiny app in the default web browser.
#' @export
#'
#' @examples
#' \donttest{
#'   run_app()
#' }
run_app <- function() {
  utils::browseURL("https://midi-pastrami.apps.math.cnrs.fr/")
}
