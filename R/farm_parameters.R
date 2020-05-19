#' Create farm_parameters object.
#'
#' @param n_farm Integer scalar. The number of farms in the simulation.
#' @param farm_area Numeric vector. Farm area (m2)
#' @param field_capacity Numeric vector. Field capacity (-)
#' @param wilting_point Numeric vector. Wilting point (-)
#' @param specific_yield Numeric vector. Specific yield (-)
#' @param runoff_coefficient Numeric vector. Runoff coefficient (-)
#' @param evaporation_coefficient Numeric vector. Evaporation coefficient (-)
#' @param canal_volume Numeric vector. Canal volume (m3)
#' @param canal_leakage_rate Numeric vector. Canal leakage rate (-)
#' @param irrigation_efficiency Numeric vector. Irrigation efficiency (-)
#'
#' @return farm_parameters
#' 
farm_parameters = function(n_farm,
                           farm_area,
                           field_capacity,
                           wilting_point,
                           specific_yield,
                           runoff_coefficient,
                           evaporation_coefficient,
                           canal_volume,
                           canal_leakage_rate,
                           irrigation_efficiency
                           ) {
    
    x = list(
        farm_area = .expand_parameter(n_farm, farm_area),
        field_capacity = .expand_parameter(n_farm, field_capacity),
        wilting_point = .expand_parameter(n_farm, wilting_point),
        specific_yield = .expand_parameter(n_farm, specific_yield),
        runoff_coefficient = .expand_parameter(n_farm, runoff_coefficient),
        evaporation_coefficient = .expand_parameter(n_farm, evaporation_coefficient),
        canal_volume = .expand_parameter(n_farm, canal_volume),
        canal_leakage_rate = .expand_parameter(n_farm, canal_leakage_rate),
        irrigation_efficiency = .expand_parameter(n_farm, canal_leakage_rate)
    )
    structure(x, class="farm_parameters")
}
