# This is where all methods (e.g. print, plot) should appear
#' @export
#' @importMethodsFrom VGAM summary
summary.compRisk <- function(object, ...) {
    summary(object$model)
}

#' @import methods
NULL
