#' Create a yield_data object.
#'
#' @param yield Data frame, xts.
#' @param time POSIX.
#'
#' @return yield_data
#' 
yield_data = function(yield, time) {
    structure(yield, "yield_data")
}
