# This is where all methods (e.g. print, plot) should appear

#' Population Time Plot
#'
#' @description \code{plot} method for objects of class \code{popTime} and \code{popTimeExposure}
#'
#' @param x an object of class \code{popTime} or \code{popTimeExposure}.
#' @param ... Ignored.
#' @param xlab,ylab,line.width,line.colour,point.size,point.colour,legend,legend.position See
#'   \code{\link[graphics]{par}}.
#' @return The methods for \code{plot} return a population time plot, stratified by exposure status
#'   in the case of \code{popTimeExposure}.
#' @import ggplot2
#' @export
#' @method plot popTime
#' @rdname popTime
plot.popTime <- function(x, ...,
                         xlab = "Follow-up time", ylab = "Population",
                         line.width = 1, line.colour = "grey80",
                         point.size = 1, point.colour = "red",
                         legend = FALSE,
                         legend.position = c("bottom", "top", "left", "right")) {
    ycoord <- yc <- `event status` <- event <- NULL

    p1 <- ggplot(x, aes(x = 0, xend = time, y = ycoord, yend = ycoord))

    p2 <- p1 +
        geom_segment(size = line.width, colour = line.colour) +
        xlab(xlab) +
        ylab(ylab) +
        theme_bw()

    if (legend) {
        legend.position <- match.arg(legend.position)
        p2 +
            geom_point(aes(x = time, y = yc, colour = `event status`),
                       data = x[event == 1], size = point.size) +
            theme(axis.text = element_text(size = 12, face = 'bold'),
                  legend.position = legend.position,
                  legend.title = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            scale_colour_manual(values = c("event" = point.colour))
    } else {
        p2 +
            geom_point(aes(x = time, y = yc),
                       data = x[event == 1], colour = point.colour,
                       size = point.size) +
            theme(axis.text = element_text(size = 12, face = 'bold'),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

    }


}

#' @param ncol Number of columns.
#' @inheritParams plot.popTime
#' @return The methods for \code{plot} return a population time plot, stratified by exposure status
#'   in the case of \code{popTimeExposure}.
#' @import ggplot2
#' @examples
#' \dontrun{
#' data(bmtccr)
#' popTimeData <- popTime(data = bmtccr, time = "ftime", exposure = "D")
#' # p is an object of class gg and ggplot
#' p <- plot(popTimeData)
#' # you can further modify the object using all ggplot2 functions
#' # here we modify the number of y-tick labels
#' p + scale_y_continuous(breaks = seq(0, max(popTimeData$data$ycoord), 10))
#' }
#' @export
#' @method plot popTimeExposure
#' @rdname popTime
plot.popTimeExposure <- function(x, ...,
                                 ncol = 1,
                                 xlab = "Follow-up time", ylab = "Population",
                                 line.width = 1, line.colour = "grey80",
                                 point.size = 1, point.colour = "red",
                                 legend = FALSE,
                                 legend.position = c("bottom", "top", "left", "right")) {


    # ds <- read.csv("data-raw/hanley/ERSPCindividualData.csv")
    # DT_ds <- as.data.table(ds)
    # DT_ds[, ScrArm:=factor(ScrArm, levels = 0:1, labels = c("No-Screening Arm","Screening Arm"))]
    #
    # object <- popTime(DT_ds, event = "DeadOfPrCa", exposure = "ScrArm")
    # roundUp(object$data[, max(ycoord)])

    # ===========================
    ycoord <- yc <- `event status` <- event <- NULL

    p1 <- ggplot(x$data, aes(x = 0, xend = time, y = ycoord, yend = ycoord))

    p2 <- p1 +
        geom_segment(size = line.width, colour = line.colour) +
        xlab(xlab) +
        ylab(ylab) +
        theme_bw() #+
        # scale_y_continuous(limits = c(0,roundUp(x$data[, max(ycoord)])))

    if (legend) {
        legend.position <- match.arg(legend.position)
        p2 +
            geom_point(aes(x = time, y = yc, colour = `event status`),
                       data = x$data[event == 1], size = point.size) +
            theme(axis.text = element_text(size = 12, face = 'bold'),
                  legend.position = legend.position,
                  legend.title = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            scale_colour_manual(values = c("event" = point.colour)) +
            facet_wrap(x$exposure, ncol = ncol)

    } else {
        p2 +
            geom_point(aes(x = time, y = yc),
                       data = x$data[event == 1], colour = point.colour,
                       size = point.size) +
            theme(axis.text = element_text(size = 12, face = 'bold'),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            facet_wrap(x$exposure, ncol = ncol)

    }


}


#' @import methods
#' @importFrom stats binomial glm integrate pnorm quantile relevel runif time update terms
NULL

#################
# S4 objects ----

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
                         eventVar = "character"),
                     contains = "vglm",
                     prototype = list(
                         originalData = data.frame(),
                         typeEvents = c(0,1),
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
setMethod("summary",
          c(object = "CompRisk"),
          function(object) callNextMethod())
