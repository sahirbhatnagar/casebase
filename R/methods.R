# S3 objects----
#' @param x Object of class \code{absRiskCB}, such as the output of
#'   \code{absoluteRisk} with single event.
#' @export
#' @rdname absoluteRisk
print.absRiskCB <- function(x, ...) {
  # Print without attributes for cleaner output
  print(x[, , drop = FALSE])
}

#################
# S4 objects----

#' @import methods
#' @importFrom stats binomial glm integrate pnorm quantile relevel runif time
#' @importFrom stats update terms
NULL

#' An S4 class to store the output of fitSmoothHazard
#'
#' This class inherits from \code{vglm-class}.
#'
#' @slot originalData Data.frame containing the original data (i.e. before
#'   case-base sampling). This is used by the \code{\link{absoluteRisk}}
#'   function.
#' @slot typeEvents Numeric factor which encodes the type of events being
#'   considered (including censoring).
#' @slot timeVar Character string giving the name of the time variable, as
#'   appearing in \code{originalData}
#' @slot eventVar Character string giving the name of the event variable, as
#'   appearing in \code{originalData}
#' @importClassesFrom VGAM vglm
CompRisk <- setClass("CompRisk",
  slots = c(
    originalData = "data.frame",
    typeEvents = "numeric",
    timeVar = "character",
    eventVar = "character"
  ),
  contains = "vglm",
  prototype = list(
    originalData = data.frame(),
    typeEvents = c(0, 1),
    timeVar = "time",
    eventVar = "event"
  )
)

#' @rdname CompRisk-class
#' @param ... Extra parameters
setGeneric("summary")
#' @export
#' @rdname CompRisk-class
#' @param object Object of class \code{CompRisk}
setMethod(
  "summary",
  c(object = "CompRisk"),
  function(object) callNextMethod()
)
