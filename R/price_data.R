#' Create a price_data object.
#'
#' @param nitrogen_price Numeric vector, xts. Nitrogen price.
#' @param phosphorous_price Numeric vector, xts. Phosphorous price.
#' @param potassium_price Numeric vector, xts. Potassium price.
#' @param diesel_price Numeric vector, xts. Diesel price.
#' @param crop_price data.frame, xts. Crop price.
#' @param time POSIX, optional.
#' @param ... Additional arguments (none)
#'
#' @return price_data
#'
price_data = function(nitrogen_price,
                      phosphorous_price,
                      potassium_price,
                      diesel_price,
                      crop_price,
                      time,
                      ...) {
    # TODO: checks
    x = data.frame(
        nitrogen_price = nitrogen_price,
        phosphorous_price = phosphorous_price,
        potassium_price = potassium_price,
        diesel_price = diesel_price,
        crop_price = crop_price
    )
    structure(x, "price_data")
}

                      
