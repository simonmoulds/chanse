#' Create an irrigation_schedule object.
#'
#' @param n_crop Integer scalar. The number of crops in the simulation.
#' @param schedule List with length equal to \code{n_crop}. Each list
#' element should be a data.frame with two columns, where the first
#' column contains the Julian days on which irrigation takes place,
#' and the second column specifies the depth of irrigation water
#' to be applied (mm).
#'
#' @return irrigation_schedule
#' 
irrigation_schedule = function(n_crop, schedule) {
    x = list(
        schedule = schedule
    )
    structure(x, class="irrigation_schedule")
}
