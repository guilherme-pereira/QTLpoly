#' Get last value of a vector
#'
#' @param x a vector from which the last value will be taken
#' @keywords internal
#' @export
last <- function(x) { tail(x, n = 1) }

#' Suppress annoying warning from varComp function
#'
#' @param w a command for fitting reduced models
#' @keywords internal
#' @export
h <- function(w) if( any( grepl( "Recycling array of length 1 in vector-array arithmetic is deprecated", w) ) ) invokeRestart( "muffleWarning" )