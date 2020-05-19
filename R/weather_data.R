#' Create a weather_data object.
#'
#' @param prec Numeric vector, xts. Daily precipitation (mm).
#' @param et0 Numeric vector, xts. Daily reference evapotranspiration (mm).
#' @param time POSIX, optional.
#'
#' @return weather_data
#'
weather_data = function(prec, et0, time) {
    ## TODO: promote to xts
    x = data.frame(
        prec = prec,
        et0 = et0
    )    
    structure(x, "weather_data")
}
