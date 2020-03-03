# This is where all methods (e.g. print, plot) should appear

#' Population Time Plot
#'
#' @description \code{plot} method for objects of class \code{popTime} and
#'   \code{popTimeExposure}
#'
#' @param x an object of class \code{popTime} or \code{popTimeExposure}.
#' @param ... Ignored.
#' @param xlab,ylab The title of the respective axis. Default: 'Follow-up time'
#'   for xlab and 'Population' for ylab
#' @param add.case.series Logical indicating if the case series should be added
#'   to the plot. Default: TRUE
#' @param add.base.series Logical indicating if the base series should be added
#'   to the plot. Default: FALSE
#' @param add.competing.event Logical indicating if the competing event should
#'   be added to the plot. Default: FALSE
#' @param casebase.theme Logical indication if the casebase theme be used. The
#'   casebase theme uses \code{\link[ggplot2]{theme_minimal}}. Default: TRUE.
#' @param ribbon.params A list containing arguments that are passed to
#'   \code{\link[ggplot2]{geom_ribbon}} which is used to plot the
#'   population-time area. These arguments will override the function defaults.
#'   For example, you can set \code{ribbon.params = list(colour = 'green')} if
#'   you want the area to be green.
#' @param case.params,base.params,competing.params A list containing arguments
#'   that are passed to \code{\link[ggplot2]{geom_point}} which is used to plot
#'   the case series, base series, competing events. These arguments will
#'   override the function defaults. For example, you can set \code{case.params
#'   = list(size = 1.5)} if you want to increase the point size for the case
#'   series points. Note: do not use this argument to change the color of the
#'   points. Doing so will result in unexpected results for the legend. See the
#'   \code{legend.params} argument, if you want to change the color of the
#'   points.
#' @param legend.params A list containing arguments that are passed to
#'   \code{\link[ggplot2]{scale_color_manual}} which is used to plot the legend.
#'   Only used if \code{legend=TRUE}. These arguments will override the function
#'   defaults. Use this argument if you want to change the color of the points.
#'   See examples for more details.
#' @param theme.params A list containing arguments that are passed to
#'   \code{\link[ggplot2]{theme}}. For example \code{theme.params =
#'   list(legend.position = 'none')}.
#' @param ratio If \code{add.base.series=TRUE}, integer, giving the ratio of the
#'   size of the base series to that of the case series. This argument is passed
#'   to the \code{\link{sampleCaseBase}} function. Default: 10.
#' @param censored.indicator If \code{add.base.series=TRUE}, a character string
#'   of length 1 indicating which value in event is the censored. This function
#'   will use relevel to set \code{censored.indicator} as the reference level.
#'   This argument is ignored if the event variable is a numeric. This argument
#'   is passed to the \code{\link{sampleCaseBase}} function.
#' @param comprisk If \code{add.base.series=TRUE}, logical indicating whether we
#'   have multiple event types and that we want to consider some of them as
#'   competing risks. This argument is passed to the
#'   \code{\link{sampleCaseBase}} function. Note: should be \code{TRUE} if your
#'   data has competing risks, even if you dont want to add competing risk
#'   points (\code{add.competing.event=FALSE}). Default: FALSE
#' @param legend Logical indicating if a legend should be added to the plot.
#'   Note that if you want to change the colors of the points, through the
#'   \code{legend.params} argument, then set \code{legend=TRUE}. If you want to
#'   change the color of the points but not have a legend, then set
#'   \code{legend=TRUE} and \code{theme.params = list(legend.position = 'none'}.
#'   Default: FALSE
#' @param legend.position Deprecated. Specify the legend.position argument
#'   instead in the \code{theme.params} argument. e.g. \code{theme.params =
#'   list(legend.position = 'bottom')}.
#' @param line.width Deprecated.
#' @param line.colour Deprecated. specify the fill argument instead in
#'   \code{ribbon.params}. e.g. \code{ribbon.params = list(fill = 'red')}.
#' @param point.size Deprecated. specify the size argument instead in the
#'   \code{case.params} or \code{base.params} or \code{competing.params}
#'   argument. e.g. \code{case.params = list(size = 1.5)}.
#' @param point.colour Deprecated. Specify the values argument instead in the
#'   \code{legend.params} argument. See examples for details.
#' @return The methods for \code{plot} return a population time plot, stratified
#'   by exposure status in the case of \code{popTimeExposure}. Note that these
#'   are \code{ggplot2} objects and can therefore be used in subsequent ggplot2
#'   type plots. See examples and vignette for details.
#' @details This function leverages the \code{ggplot2} package to build
#'   population time plots. It builds the plot by adding layers, starting with a
#'   layer for the area representing the population time. It then sequentially
#'   adds points to the plots to show the casebase sampling mechanism. This
#'   function gives user the flexibility to add any combination of the
#'   case.series, base.series and competing events. The case series and
#'   competing events are sampled at random vertically on the plot in order to
#'   visualise the incidence density using the \code{\link{popTime}} function.
#'   That is, imagine we draw a vertical line at a specific event time. We then
#'   plot the point at a randomly sampled y-coordinate along this vertical line.
#'   This is done to avoid having all points along the upper edge of the plot
#'   (because the subjects with the least amount of observation time are plotted
#'   at the top of the y-axis). By randomly distributing them, we can get a
#'   better sense of the inicidence density. The base series is sampled
#'   horizontally on the plot using the \code{\link{sampleCaseBase}} function.
#' @import ggplot2
#' @seealso
#'   \link[ggplot2]{geom_point},\link[ggplot2]{geom_ribbon},\link[ggplot2]{theme},
#'   \link[ggplot2]{scale_colour_manual}, \link{sampleCaseBase}
#' @examples
#' # change color of points, but don't produce a legend
#' library(ggplot2)
#' data("bmtcrr")
#' popTimeData <- popTime(data = bmtcrr, time = "ftime", event = "Status")
#' cols <- c("Case series" = "black", "Competing event" = "#009E73", "Base series" = "#0072B2")
#' plot(popTimeData,
#'      casebase.theme = TRUE,
#'      add.case.series = TRUE,
#'      ratio = 1,
#'      add.base.series = TRUE,
#'      legend = TRUE,
#'      comprisk = TRUE,
#'      add.competing.event = FALSE,
#'      legend.params = list(name = element_blank(),
#'                      breaks = c("Case series", "Competing event", "Base series"),
#'                      values = cols),
#'      theme.params = list(legend.position = "none"))
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
                         ribbon.params = list(),
                         case.params = list(),
                         base.params = list(),
                         competing.params = list(),
                         legend.params = list(),
                         theme.params = list(),
                         ratio = 10,
                         censored.indicator,
                         comprisk = FALSE,
                         legend = FALSE,
                         legend.position,
                         line.width,
                         line.colour,
                         point.size,
                         point.colour) {

    ycoord <- yc <- `event status` <- event <- comprisk.event <- NULL

    if (!missing(line.colour)) {
        warning("line.colour argument deprecated. specify the fill argument instead
                in the ribbon.params. e.g. ribbon.params = list(fill = 'red').")
    }

    if (!missing(line.width)) {
        warning("line.width argument deprecated.")
    }

    if (!missing(point.size)) {
        warning("point.size argument deprecated. specify the size argument instead
                in the case.params or base.params or competing.params argument.
                e.g. case.params = list(size = 1.5).")
    }

    if (!missing(point.colour)) {
        warning("point.colour argument deprecated. specify the values argument instead
                in the legend.params argument. see examples for details.")
    }

    if (!missing(legend.position)) {
        warning("legend.position argument deprecated. specify the legend.position argument instead
                in the theme.params argument.
                e.g. theme.params = list(legend.position = 'bottom').")
    }

    p2 <- list(

        # Add poptime area --------------------------------------------------------
        do.call("geom_ribbon", utils::modifyList(
            list(data = x,
                 mapping = aes(x = time, ymin = 0, ymax = ycoord),
                 fill = "grey80"),
            ribbon.params)
        ),

        # Add case series ---------------------------------------------------------
        if (add.case.series)
            do.call("geom_point", utils::modifyList(
                list(data = x[event == 1],
                     mapping = aes(x = time, y = yc, colour = "Case series"),
                     show.legend = legend),
                case.params)
            ),


        # Add base series ---------------------------------------------------------
        if (add.base.series) {
            basedata <- sampleCaseBase(data = x, time = "time", event = "event", ratio = ratio,
                                       comprisk = comprisk,
                                       censored.indicator = censored.indicator)
            # browser()
            do.call("geom_point", utils::modifyList(
                list(data = basedata[event == 0],
                     mapping = aes(x = time, y = ycoord, colour = "Base series"),
                     show.legend = legend),
                base.params)
            )
        },


        # Add competing event -----------------------------------------------------
        if (add.competing.event) {

            newX <- data.table::copy(x)
            newX[, `:=`(original.time = NULL,
                        original.event = NULL,
                        `event status` = NULL,
                        ycoord = NULL,
                        yc = NULL,
                        n_available = NULL)]
            newX[event == 0, comprisk.event := 0]
            newX[event == 1, comprisk.event := 2]
            newX[event == 2, comprisk.event := 1]
            newX[, event := NULL]

            compdata <- popTime(data = newX, time = "time", event = "comprisk.event")

            do.call("geom_point", utils::modifyList(
                list(data = compdata[event == 1],
                     mapping = aes(x = time, y = yc, colour = "Competing event"),
                     show.legend = legend),
                competing.params)
            )
        },

        # Add legend --------------------------------------------------------------
        if (legend) {

            cols <- c("Case series" = "#D55E00", "Competing event" = "#009E73", "Base series" = "#0072B2")

            do.call("scale_colour_manual", utils::modifyList(
                list(name = element_blank(),
                     breaks = c("Case series", "Competing event", "Base series"),
                     values = cols), legend.params)
            )
        },

        # Use casebase theme or not -----------------------------------------------
        if (casebase.theme)
            theme_minimal(),


        # add theme stuff  --------------------------------------------------------
        do.call("theme", theme.params)
    )


    # cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # plot(seq_along(cbPalette),col = cbPalette, pch = 19, cex = 2.5)

    p1 <- ggplot()

    p1 + p2 + xlab(xlab) + ylab(ylab)
}



#' @param ncol Deprecated. Use \code{facet.params} instead.
#' @param facet.params A list containing arguments that are passed to
#'   \code{\link[ggplot2]{facet_wrap}} which is used to create facet plots. Only
#'   used if plotting exposure stratified population time plots. These arguments
#'   will override the function defaults.
#' @return The methods for \code{plot} return a population time plot, stratified
#'   by exposure status in the case of \code{popTimeExposure}.
#' @import ggplot2
#' @examples
#' \dontrun{
#' data("bmtcrr")
#' popTimeData <- popTime(data = bmtcrr, time = "ftime", exposure = "D")
#' # p is an object of class gg and ggplot
#' p <- plot(popTimeData)
#' # you can further modify the object using all ggplot2 functions
#' # here we modify the number of y-tick labels
#' p + scale_y_continuous(breaks = seq(0, max(popTimeData$ycoord), 10))
#' }
#' @export
#' @method plot popTimeExposure
#' @rdname popTime
plot.popTimeExposure <- function(x, ...,
                                 xlab = "Follow-up time",
                                 ylab = "Population",
                                 add.case.series = TRUE,
                                 add.base.series = FALSE,
                                 add.competing.event = FALSE,
                                 casebase.theme = TRUE,
                                 ribbon.params = list(),
                                 case.params = list(),
                                 base.params = list(),
                                 competing.params = list(),
                                 legend.params = list(),
                                 theme.params = list(),
                                 facet.params = list(),
                                 ratio = 10,
                                 censored.indicator,
                                 comprisk = FALSE,
                                 legend = FALSE,
                                 ncol,
                                 legend.position,
                                 line.width,
                                 line.colour,
                                 point.size,
                                 point.colour) {

    # ===========================
    ycoord <- yc <- `event status` <- event <- comprisk.event <- NULL

    exposure_variable <- attr(x, "exposure")

    if (!missing(line.colour)) {
        warning("line.colour argument deprecated. specify the fill argument instead
                in the ribbon.params. e.g. ribbon.params = list(fill = 'red').")
    }

    if (!missing(line.width)) {
        warning("line.width argument deprecated.")
    }

    if (!missing(point.size)) {
        warning("point.size argument deprecated. specify the size argument instead
                in the case.params or base.params or competing.params argument.
                e.g. case.params = list(size = 1.5).")
    }

    if (!missing(ncol)) {
        warning("ncol argument deprecated. specify the ncol argument instead
                in the facet.params argument.
                e.g. facet.params = list(ncol = 1).")
    }

    if (!missing(point.colour)) {
        warning("point.colour argument deprecated. specify the values argument instead
                in the legend.params argument. see examples for details.")
    }

    if (!missing(legend.position)) {
        warning("legend.position argument deprecated. specify the legend.position argument instead
                in the theme.params argument.
                e.g. theme.params = list(legend.position = 'bottom').")
    }

    p2 <- list(

        # Add poptime area --------------------------------------------------------
        do.call("geom_ribbon", utils::modifyList(
            list(data = x,
                 mapping = aes(x = time, ymin = 0, ymax = ycoord),
                 fill = "grey80"),
            ribbon.params)
        ),

        # Add case series ---------------------------------------------------------
        if (add.case.series)
            do.call("geom_point", utils::modifyList(
                list(data = x[event == 1],
                     mapping = aes(x = time, y = yc, colour = "Case series"),
                     show.legend = legend),
                case.params)
            ),


        # Add base series ---------------------------------------------------------
        if (add.base.series) {
            basedata <- sampleCaseBase(data = x, time = "time", event = "event", ratio = ratio,
                                       comprisk = comprisk,
                                       censored.indicator = censored.indicator)
            # browser()
            do.call("geom_point", utils::modifyList(
                list(data = basedata[event == 0],
                     mapping = aes(x = time, y = ycoord, colour = "Base series"),
                     show.legend = legend),
                base.params)
            )
        },


        # Add competing event -----------------------------------------------------
        if (add.competing.event) {

            newX <- data.table::copy(x)
            newX[, `:=`(original.time = NULL,
                        original.event = NULL,
                        `event status` = NULL,
                        ycoord = NULL,
                        yc = NULL,
                        n_available = NULL)]
            newX[event == 0, comprisk.event := 0]
            newX[event == 1, comprisk.event := 2]
            newX[event == 2, comprisk.event := 1]
            newX[, event := NULL]

            compdata <- popTime(data = newX, time = "time", event = "comprisk.event",
                                exposure = exposure_variable)

            do.call("geom_point", utils::modifyList(
                list(data = compdata[event == 1],
                     mapping = aes(x = time, y = yc, colour = "Competing event"),
                     show.legend = legend),
                competing.params)
            )
        },

        # Add legend --------------------------------------------------------------
        if (legend) {

            cols <- c("Case series" = "#D55E00", "Competing event" = "#009E73", "Base series" = "#0072B2")

            do.call("scale_colour_manual", utils::modifyList(
                list(name = element_blank(),
                     breaks = c("Case series", "Competing event", "Base series"),
                     values = cols), legend.params)
            )
        },

        do.call("facet_wrap", utils::modifyList(
            list(facets = exposure_variable, ncol = 1),
            facet.params
        )),

        # Use casebase theme or not -----------------------------------------------
        if (casebase.theme)
            theme_minimal(),


        # add theme stuff  --------------------------------------------------------
        do.call("theme", theme.params)
    )


    # cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # plot(seq_along(cbPalette),col = cbPalette, pch = 19, cex = 2.5)

    p1 <- ggplot()

    p1 + p2 + xlab(xlab) + ylab(ylab)

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
