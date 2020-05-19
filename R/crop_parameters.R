#' Create a crop_parameters object.
#'
#' @param n_crop Integer scalar. The number of crops in the simulation.
#' @param L_ini Numeric vector. The relative length of the inital crop
#' growth phase.
#' @param L_dev Numeric vector. The relative length of the development
#' crop growth phase.
#' @param L_mid Numeric vector.
#' @param L_late Numeric vector.
#' @param Kc_ini Numeric vector.
#' @param Kc_mid Numeric vector.
#' @param Kc_end Numeric vector.
#' @param rd Numeric vector.
#'
#' @return crop_parameters
#' 
crop_parameters = function(n_crop,
                           L_ini,
                           L_dev,
                           L_mid,
                           L_late,
                           Kc_ini,
                           Kc_mid,
                           Kc_end,
                           rd) {
    x = list(
        L_ini = .expand_parameter(n_crop, L_ini),
        L_dev = .expand_parameter(n_crop, L_dev),
        L_mid = .expand_parameter(n_crop, L_mid),
        L_late = .expand_parameter(n_crop, L_late),
        Kc_ini = .expand_parameter(n_crop, Kc_ini),
        Kc_mid = .expand_parameter(n_crop, Kc_mid),
        Kc_end = .expand_parameter(n_crop, Kc_end),
        rd = .expand_parameter(n_crop, rd)
    )
    structure(x, class="crop_parameters")
}                           
