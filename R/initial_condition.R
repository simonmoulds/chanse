#' Create initial_condition object.
#'
#' @param n_farm Integer scalar. The number of farms in the model simulation.
#' @param root_zone_depletion Numeric vector. Initial root zone depletion.
#' @param groundwater_depth Numeric vector. Initial depth to groundwater.
#' @param initial_savings Numeric vector. Initial farmer savings.
#' @param initial_category Integer vector. Initial category.
#'
#' @return initial_condition
#' 
initial_condition = function(n_farm,
                             root_zone_depletion,
                             groundwater_depth,
                             initial_savings,
                             initial_category
                             ) {
    
    x = list(
        root_zone_depletion = .expand_parameter(n_farm, root_zone_depletion),
        groundwater_depth = .expand_parameter(n_farm, groundwater_depth),
        savings = .expand_parameter(n_farm, initial_savings),
        category = .expand_parameter(n_farm, initial_category)
    )
    structure(x, class="initial_condition")
}
