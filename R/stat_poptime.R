#' currently not being used
#' p <- ggplot(DT, aes(x=0, xend=time, y=ycoord, yend=ycoord))
#'
#' p + theme_classic() +
#'     geom_segment(size=3, colour="grey90") +
#'     geom_point(aes(x=time, y=yc), data = DT[event==1], size=1, colour="red") +
#'     theme(axis.text=element_text(size=12, face='bold'))
#'
#' geom_xspline <- function(mapping = NULL, data = NULL, stat = "xspline",
#'                          position = "identity", show.legend = NA,
#'                          inherit.aes = TRUE, na.rm = TRUE,
#'                          spline_shape=-0.25, open=TRUE, rep_ends=TRUE, ...) {
#'     layer(
#'         geom = GeomXspline,
#'         mapping = mapping,
#'         data = data,
#'         stat = stat,
#'         position = position,
#'         show.legend = show.legend,
#'         inherit.aes = inherit.aes,
#'         params = list(spline_shape=spline_shape,
#'                       open=open,
#'                       rep_ends=rep_ends,
#'                       ...)
#'     )
#' }
#'
#' GeomXspline <- ggproto("GeomXspline", GeomSegment,
#'                        required_aes = c("x", "xend", "y", "yend"),
#'                        default_aes = aes(colour = "grey90", size = 1,
#'                                          linetype = 1)
#' )
#'
#'
#'
#' #' Empirical Receiver Operating Characteristic Curve
#' #'
#' #' Display the empirical ROC curve. Useful for characterizing the classification
#' #' accuracy of continuous measurements for predicting binary states
#' #'
#' #' @section Aesthetics:
#' #' \code{geom_roc} understands the following aesthetics (required aesthetics
#' #' are in bold):
#' #' \itemize{
#' #'   \item \strong{\code{x}} The FPF estimate. This is automatically mapped by \link{stat_roc}
#' #'   \item \strong{\code{y}} The TPF estimate. This is automatically mapped by \link{stat_roc}
#' #'   smallest level in sort order is assumed to be 0, with a warning
#' #'   \item \code{alpha}
#' #'   \item \code{color}
#' #'   \item \code{fill}
#' #'   \item \code{linetype}
#' #'   \item \code{size}
#' #' }
#' #'
#' #' @param stat Use to override the default connection between
#' #'   \code{geom_roc} and \code{stat_roc}.
#' #' @seealso See \code{\link{geom_rocci}} for
#' #'   displaying rectangular confidence regions for the empirical ROC curve, \code{\link{style_roc}} for
#' #'   adding guidelines and labels, and \code{\link{direct_label}} for adding direct labels to the
#' #'   curves. Also \link{export_interactive_roc} for creating interactive ROC curve plots for use in a web browser.
#' #' @inheritParams ggplot2::geom_point
#' #' @export
#' #'
#'
#' geom_roc <- function(mapping = NULL, data = NULL, stat = "roc", n.cuts = 10, arrow = NULL,
#'                      lineend = "butt", linejoin = "round", linemitre = 1,
#'                      alpha.line = 1, alpha.point = 1,
#'                      size.point = .5, labels = TRUE, labelsize = 3.88, labelround = 1,
#'                      na.rm = TRUE, position = "identity", show.legend = NA, inherit.aes = TRUE, ...) {
#'     layer(
#'         geom = GeomRoc, mapping = mapping, data = data, stat = stat,
#'         position = position, show.legend = show.legend, inherit.aes = inherit.aes,
#'         params = list(na.rm = na.rm, n.cuts = n.cuts, arrow = arrow,
#'                       lineend = lineend, linejoin = linejoin, linemitre = linemitre,
#'                       alpha.line = alpha.line, alpha.point = alpha.point,
#'                       size.point = size.point, labels = labels, labelsize = labelsize, labelround = labelround,
#'                       ...)
#'     )
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' StatLm <- ggproto("StatLm", Stat,
#'                   required_aes = c("x", "y"),
#'
#'                   compute_group = function(data, scales, params, n = 100, formula = y ~ x) {
#'                       rng <- range(data$x, na.rm = TRUE)
#'                       grid <- data.frame(x = seq(rng[1], rng[2], length = n))
#'
#'                       mod <- lm(formula, data = data)
#'                       grid$y <- predict(mod, newdata = grid)
#'
#'                       grid
#'                   }
#' )
#'
#' stat_lm <- function(mapping = NULL, data = NULL, geom = "line",
#'                     position = "identity", na.rm = FALSE, show.legend = NA,
#'                     inherit.aes = TRUE, n = 50, formula = y ~ x,
#'                     ...) {
#'     layer(
#'         stat = StatLm, data = data, mapping = mapping, geom = geom,
#'         position = position, show.legend = show.legend, inherit.aes = inherit.aes,
#'         params = list(n = n, formula = formula, na.rm = na.rm, ...)
#'     )
#' }
#'
#'
