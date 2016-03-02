# This is where all methods (e.g. print, plot) should appear



#' Population Time Plot
#'
#' \code{plot} method for objects of class \code{popTime}
#'
#' @param object an object of class \code{popTime}. See \code{\link{popTime}}
#'   for details
#' @return a population time plot
#' @import ggplot2
#' @export
plot.popTime <- function(object,
                         xlab = "Follow-up years", ylab = "Population",
                         line.width = 1, line.colour = "grey80",
                         point.size = 1, point.colour = "red",
                         legend = FALSE,
                         legend.position = c("bottom", "top", "left", "right")) {

    p1 <- ggplot(object, aes(x=0, xend=time, y=ycoord, yend=ycoord))

    p2 <- p1 +
        geom_segment(size = line.width, colour = line.colour) +
        xlab(xlab) +
        ylab(ylab) +
        theme_bw()

    if (legend) {
        legend.position <- match.arg(legend.position)
        p2 +
            geom_point(aes(x=time, y=yc, colour = `event status`),
                       data = object[event==1], size = point.size) +
            theme(axis.text=element_text(size=12, face='bold'),
                  legend.position = legend.position,
                  legend.title=element_blank()) +
            scale_colour_manual(values = c("event" = point.colour))
    } else {
        p2 +
            geom_point(aes(x=time, y=yc),
                       data = object[event==1], colour = point.colour,
                       size = point.size) +
            theme(axis.text=element_text(size=12, face='bold'))

    }


}





#' Population Time Plot Stratified by Exposure Status
#'
#' \code{plot} method for objects of class \code{popTimeExposure}
#'
#' @param object an object of class \code{popTimeExposure}. See
#'   \code{\link{popTime}} for details
#' @return a population time plot stratified by exposure status
#' @import ggplot2
#' @export
plot.popTimeExposure <- function(object,
                                 ncol = 2,
                                 xlab = "Follow-up years", ylab = "Population",
                                 line.width = 1, line.colour = "grey80",
                                 point.size = 1, point.colour = "red",
                                 legend = FALSE,
                                 legend.position = c("bottom", "top", "left", "right")) {

    p1 <- ggplot(object$data, aes(x=0, xend=time, y=ycoord, yend=ycoord))

    p2 <- p1 +
        geom_segment(size = line.width, colour = line.colour) +
        xlab(xlab) +
        ylab(ylab) +
        theme_bw()

    if (legend) {
        legend.position <- match.arg(legend.position)
        p2 +
            geom_point(aes(x=time, y=yc, colour = `event status`),
                       data = object$data[event==1], size = point.size) +
            theme(axis.text=element_text(size=12, face='bold'),
                  legend.position = legend.position,
                  legend.title=element_blank()) +
            scale_colour_manual(values = c("event" = point.colour)) +
            facet_wrap(object$exposure, ncol = ncol)
    } else {
        p2 +
            geom_point(aes(x=time, y=yc),
                       data = object$data[event==1], colour = point.colour,
                       size = point.size) +
            theme(axis.text=element_text(size=12, face='bold')) +
            facet_wrap(object$exposure, ncol = ncol)

    }


}
