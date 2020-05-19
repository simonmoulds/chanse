#' @export
#' @noRd
.expand_parameter = function(n, param) {
    if (length(param) == n) {
        return(param)
    } else if (length(param) == 1) {
        return(rep(param, n))
    } else {
        stop()
    }
}
