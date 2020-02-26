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
                         xlab = "Follow-up time",
                         ylab = "Population",
                         add.case.series = TRUE,
                         add.base.series = FALSE,
                         add.competing.event = FALSE,
                         casebase.theme = TRUE,
                         area.colour = "grey80",
                         ribbon.params = list(),
                         case.params = list(),
                         base.params = list(),
                         competing.params = list(),
                         legend.params = list(),
                         ratio = 10,
                         comprisk = FALSE,
                         censored.indicator,
                         line.width,
                         line.colour,
                         point.size = 1,
                         point.colour = "red",
                         legend = FALSE,
                         legend.position = c("bottom", "top", "left", "right")) {

    ycoord <- yc <- `event status` <- event <- NULL

    if (!missing(line.colour)) {
        warning("line.colour argument deprecated. use area.colour argument instead.
                The parameter area.colour is set equal to line.colour.")
        area.colour <- line.colour
    }

    if (!missing(line.width)) {
        warning("line.width argument deprecated.")
    }
browser()
    params <- list(...)
    case.params <- modifyList(params, case.params)
    base.params  <- modifyList(params, base.params)
    competing.params  <- modifyList(params, competing.params)
    ribbon.params  <- modifyList(params, ribbon.params)

    p2 <- list(

        # Add poptime area --------------------------------------------------------
        do.call("geom_ribbon", modifyList(
            list(data = x, mapping = aes(x = time, ymin = 0, ymax = ycoord), fill = "grey80"),
            ribbon.params)
        ),

        # Add case series ---------------------------------------------------------
        if (add.case.series)
            do.call("geom_point", modifyList(
                list(data = x[event == 1], mapping = aes(x = time, y = yc, colour = "Case series")),
                case.params)
            ),


        # Add base series ---------------------------------------------------------
        if (add.base.series) {
            basedata <- sampleCaseBase(data = x, time = "time", event = "event", ratio = ratio,
                                       comprisk = comprisk, censored.indicator = censored.indicator)
            do.call("geom_point", modifyList(
                list(data = basedata[event == 0], mapping = aes(x = time, y = ycoord, colour = "Base series")),
                base.params)
            )
        },


        # Add legend --------------------------------------------------------------
        if (legend) {
            cols <- c("Case series" = "#D55E00", "Competing event" = "#009E73", "Base series" = "#0072B2")

            # cols <- c("Case series" = "#D55E00", "Base series" = "#0072B2")


            do.call("scale_colour_manual", modifyList(
                list(name = element_blank(),
                     breaks = c("Case series", "Competing evemt", "Base series"),
                     # breaks = c("Case series", "Base series"),
                     values = cols), legend.params)
            )
        },


        # Use casebase theme or not -----------------------------------------------
        if (casebase.theme)
            theme_minimal()
    )


    # cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # plot(seq_along(cbPalette),col = cbPalette, pch = 19, cex = 2.5)

browser()

    # p1 <- ggplot(x, aes(x = 0, xend = time, y = ycoord, yend = ycoord))
    p1 <- ggplot()

    p1 + p2 + xlab(xlab) + ylab(ylab)
    # browser()
    # #+
    #     # scale_colour_manual(name = element_blank(),
    #                        # breaks = c("Case series", "Base series"),
    #                        # values = cols) #+
    #     # theme(axis.text = element_text(size = 12, face = 'bold'),
    #           # legend.position = "bottom",
    #           # legend.title = element_blank())
    #
    #
    # if (legend) {
    #     legend.position <- match.arg(legend.position)
    #     p2 +
    #         geom_point(aes(x = time, y = yc, colour = `event status`),
    #                    data = x[event == 1], size = point.size) +
    #         theme(axis.text = element_text(size = 12, face = 'bold'),
    #               legend.position = legend.position,
    #               legend.title = element_blank())
    # } else {
    #     p2 +
    #         geom_point(aes(x = time, y = yc),
    #                    data = x[event == 1], colour = point.colour,
    #                    size = point.size) +
    #         theme(axis.text = element_text(size = 12, face = 'bold'))#,
    #               # panel.grid.major = element_blank(),
    #               # panel.grid.minor = element_blank())
    #
    # }


}



#' @param ncol Number of columns.
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
