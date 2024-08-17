#' Population Time Plot
#'
#' @description \code{plot} method for objects of class \code{popTime}
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
#'   casebase theme uses [ggplot2::theme_minimal()]. Default: TRUE.
#' @param ribbon.params A list containing arguments that are passed to
#'   [ggplot2::geom_ribbon()] which is used to plot the
#'   population-time area. These arguments will override the function defaults.
#'   For example, you can set \code{ribbon.params = list(colour = 'green')} if
#'   you want the area to be green.
#' @param case.params,base.params,competing.params A list containing arguments
#'   that are passed to [ggplot2::geom_point()] which is used to plot
#'   the case series, base series, competing events. These arguments will
#'   override the function defaults. For example, you can set \code{case.params
#'   = list(size = 1.5)} if you want to increase the point size for the case
#'   series points. Note: do not use this argument to change the color of the
#'   points. Doing so will result in unexpected results for the legend. See the
#'   \code{color.params} and \code{fill.params} arguments, if you want to change
#'   the color of the points.
#' @param color.params A list containing arguments that are passed to
#'   [ggplot2::scale_color_manual()] which is used to plot the legend.
#'   Only used if \code{legend=TRUE}. These arguments will override the function
#'   defaults. Use this argument if you want to change the color of the points.
#'   See examples for more details.
#' @param fill.params A list containing arguments that are passed to
#'   [ggplot2::scale_fill_manual()] which is used to plot the legend.
#'   Only used if \code{legend=TRUE}. These arguments will override the function
#'   defaults. Use this argument if you want to change the color of the points.
#'   See examples for more details.
#' @param theme.params A list containing arguments that are passed to
#'   [ggplot2::theme()]. For example \code{theme.params =
#'   list(legend.position = 'none')}.
#' @param facet.params A list containing arguments that are passed to
#'   [ggplot2::facet_wrap()] which is used to create facet plots. Only
#'   used if plotting exposure stratified population time plots. These arguments
#'   will override the function defaults.
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
#'   data has competing risks, even if you don't want to add competing risk
#'   points (\code{add.competing.event=FALSE}). Default: FALSE
#' @param legend Logical indicating if a legend should be added to the plot.
#'   Note that if you want to change the colors of the points, through the
#'   \code{color.params} and \code{fill.params} arguments, then set
#'   \code{legend=TRUE}. If you want to change the color of the points but not
#'   have a legend, then set \code{legend=TRUE} and \code{theme.params =
#'   list(legend.position = 'none'}. Default: FALSE
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
#'   \code{color.params} and \code{fill.params} argument. See examples for
#'   details.
#' @param ncol Deprecated. Use \code{facet.params} instead.
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
#'   better sense of the incidence density. The base series is sampled
#'   horizontally on the plot using the \code{\link{sampleCaseBase}} function.
#' @importFrom data.table := copy
#' @importFrom ggplot2 ggplot geom_point scale_fill_manual geom_ribbon aes
#' @importFrom ggplot2 scale_colour_manual element_blank facet_wrap
#' @importFrom ggplot2 theme xlab ylab
#' @seealso
#' [ggplot2::geom_point()],[ggplot2::geom_ribbon()],[ggplot2::theme()],
#' [ggplot2::scale_colour_manual()], [ggplot2::scale_fill_manual()],
#' \link{sampleCaseBase}
#' @examples
#' # change color of points
#' library(ggplot2)
#' data("bmtcrr")
#' popTimeData <- popTime(data = bmtcrr, time = "ftime", event = "Status")
#' fill_cols <- c("Case series" = "black", "Competing event" = "#009E73",
#'                "Base series" = "#0072B2")
#' color_cols <- c("Case series" = "black", "Competing event" = "black",
#'                 "Base series" = "black")
#'
#' plot(popTimeData,
#'   add.case.series = TRUE,
#'   add.base.series = TRUE,
#'   add.competing.event = FALSE,
#'   legend = TRUE,
#'   comprisk = TRUE,
#'   fill.params = list(
#'     name = element_blank(),
#'     breaks = c("Case series", "Competing event", "Base series"),
#'     values = fill_cols
#'   ),
#'   color.params = list(
#'     name = element_blank(),
#'     breaks = c("Case series", "Competing event", "Base series"),
#'     values = color_cols
#'   )
#' )
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
                         color.params = list(),
                         fill.params = list(),
                         theme.params = list(),
                         facet.params = list(),
                         ratio = 1,
                         censored.indicator,
                         comprisk = FALSE,
                         legend = TRUE,
                         ncol, # Deprecated
                         legend.position, # Deprecated
                         line.width, # Deprecated
                         line.colour, # Deprecated
                         point.size, # Deprecated
                         point.colour) { # Deprecated

    # To prevent "no visble binding for global variable"
    ycoord <- yc <- event <- comprisk.event <- NULL

    # Okabe Ito colors
    # fill_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    #                "#999999")
    # colorspace::darken(fill_cols, 0.3)
    # color_cols <- c("#9D6C06", "#077DAA", "#026D4E", "#A39A09", "#044F7E",
    #                 "#696969")

    # cbbPalette
    # fill_cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    #                "#0072B2", "#D55E00", "#CC79A7")
    # colorspace::qualitative_hcl(n = 3, palette = "Dark3") %>% put
    fill_cols <- c("#E16A86", "#50A315", "#009ADE")
    # colorspace::darken(fill_cols, 0.3) %>% dput
    # color_cols <- c("#696969", "#9D6C06", "#077DAA", "#026D4E", "#A39A09",
    #                 "#044F7E", "#954000", "#984B77")
    color_cols <- c("#AB3A59", "#347004", "#026A9A")
    # par(mfrow=c(1,2))
    # plot(seq_along(color_cols),col = color_cols, pch = 19, cex = 2.5)
    # plot(seq_along(fill_cols),col = fill_cols, pch = 19, cex = 2.5)

    if (!missing(line.colour)) {
        warning(paste("line.colour argument deprecated. specify the fill",
                      "argument instead in the ribbon.params.",
                      "e.g. ribbon.params = list(fill = 'red')."))
    }

    if (!missing(line.width)) {
        warning("line.width argument deprecated.")
    }

    if (!missing(point.size)) {
        warning(paste("point.size argument deprecated. specify the size",
                      "argument instead in the case.params or base.params or",
                      "competing.params argument.",
                      "e.g. case.params = list(size = 1.5)."))
    }

    if (!missing(point.colour)) {
        warning(paste("point.colour argument deprecated. specify the values",
                      "argument instead in the fill.params and color.params",
                      "arguments. see examples for details."))
    }

    if (!missing(legend.position)) {
        warning(paste("legend.position argument deprecated. specify the",
                      "legend.position argument instead in the theme.params",
                      "argument.",
                      "e.g. theme.params = list(legend.position = 'bottom')."))
    }

    if (!missing(ncol)) {
        warning(paste("ncol argument deprecated. specify the ncol argument",
                      "instead in the facet.params argument.",
                      "e.g. facet.params = list(ncol = 1)."))
    }

    if (length(fill.params) > 0 & length(color.params) == 0) {
        warning(paste("fill.params has been specified by the user but",
                      "color.params has not. Setting color.params to be equal",
                      "to fill.params."))
        color.params <- fill.params
    }

    if (length(color.params) > 0 & length(fill.params) == 0) {
        warning(paste("color.params has been specified by the user but",
                      "fill.params has not. Setting fill.params to be equal",
                      "to color.params."))
        fill.params <- color.params
    }

    exposure_variable <- attr(x, "exposure")

    p2 <- list(

        # Add poptime area -----------------------------------------------------
        do.call("geom_ribbon", utils::modifyList(
            list(
                data = x,
                mapping = aes(x = time, ymin = 0, ymax = ycoord),
                fill = "grey80",
                alpha = 0.5
            ),
            ribbon.params
        )),

        # Add case series ------------------------------------------------------
        if (add.case.series) {
            do.call("geom_point", utils::modifyList(
                list(
                    data = x[event == 1],
                    mapping = aes(x = time, y = yc, color = "Case series",
                                  fill = "Case series"),
                    show.legend = legend,
                    size = 1.5,
                    alpha = 0.5,
                    shape = 21
                ),
                case.params
            ))
        },


        # Add base series ------------------------------------------------------
        if (add.base.series) {
            basedata <- sampleCaseBase(
                data = x, time = "time", event = "event", ratio = ratio,
                comprisk = comprisk,
                censored.indicator = censored.indicator
            )

            do.call("geom_point", utils::modifyList(
                list(
                    data = basedata[event == 0],
                    mapping = aes(x = time, y = ycoord, colour = "Base series",
                                  fill = "Base series"),
                    show.legend = legend,
                    size = 1.5,
                    alpha = 0.5,
                    shape = 21
                ),
                base.params
            ))
        },


        # Add competing event --------------------------------------------------
        if (add.competing.event) {
            newX <- data.table::copy(x)
            newX[, `:=`(
                original.time = NULL,
                original.event = NULL,
                `event status` = NULL,
                ycoord = NULL,
                yc = NULL,
                n_available = NULL
            )]
            newX[event == 0, comprisk.event := 0]
            newX[event == 1, comprisk.event := 2]
            newX[event == 2, comprisk.event := 1]
            newX[, event := NULL]

            if (is.null(exposure_variable)) {
                compdata <- popTime(data = newX, time = "time",
                                    event = "comprisk.event")
            } else {
                compdata <- popTime(
                    data = newX, time = "time", event = "comprisk.event",
                    exposure = exposure_variable
                )
            }

            do.call("geom_point", utils::modifyList(
                list(
                    data = compdata[event == 1],
                    mapping = aes(x = time, y = yc, colour = "Competing event",
                                  fill = "Competing event"),
                    show.legend = legend,
                    size = 1.5,
                    alpha = 0.5,
                    shape = 21
                ),
                competing.params
            ))
        },

        # Add legend -----------------------------------------------------------
        if (legend) {

            # cols <- c("Case series" = color_cols[7],
            #           "Competing event" = color_cols[4],
            #           "Base series" = color_cols[6])

            cols <- c(
                "Case series" = color_cols[1],
                "Competing event" = color_cols[3],
                "Base series" = color_cols[2]
            )

            do.call("scale_colour_manual", utils::modifyList(
                list(
                    name = element_blank(),
                    breaks = c("Case series", "Competing event", "Base series"),
                    values = cols
                ), color.params
            ))
        },

        if (legend) {

            # cols <- c("Case series" = fill_cols[7],
            #           "Competing event" = fill_cols[4],
            #           "Base series" = fill_cols[6])

            cols <- c(
                "Case series" = fill_cols[1],
                "Competing event" = fill_cols[3],
                "Base series" = fill_cols[2]
            )

            do.call("scale_fill_manual", utils::modifyList(
                list(
                    name = element_blank(),
                    breaks = c("Case series", "Competing event", "Base series"),
                    values = cols
                ), fill.params
            ))
        },


        # Exposure stratified --------------------------------------------------
        if (!is.null(exposure_variable)) {
            do.call("facet_wrap", utils::modifyList(
                list(facets = exposure_variable, ncol = 1),
                facet.params
            ))
        },

        # Use casebase theme or not --------------------------------------------
        if (casebase.theme) {
            if (!is.null(exposure_variable)) {
                theme_cb() + panelBorder()
            } else {
                theme_cb()
            }
        },

        # add theme stuff  -----------------------------------------------------
        do.call("theme", utils::modifyList(
            list(legend.position = "bottom"),
            theme.params
        ))
    )

    p1 <- ggplot()

    p1 + p2 + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
}



#' @title Plot Hazards and Hazard Ratios
#'
#' @description Plot method for objects returned by the \code{fitSmoothHazard}
#'   function. Current plot types are hazard function and hazard ratio. The
#'   \code{visreg} package must be installed for \code{type="hazard"}. This
#'   function accounts for the possible time-varying exposure effects.
#' @param x Fitted object of class `glm`, `gam`, `cv.glmnet` or `gbm`. This is
#'   the result from the [casebase::fitSmoothHazard()] function.
#' @param ... further arguments passed to `plot`. Only used if \code{type="hr"}.
#'   Any of `lwd`,`lty`,`col`,`pch`,`cex` will be applied to the hazard ratio
#'   line, or point (if only one time point is supplied to `newdata`).
#' @param type plot type. Choose one of either \code{"hazard"} for hazard
#'   function or \code{"hr"} for hazard ratio.  Default: \code{type = "hazard"}.
#' @param hazard.params Named list of arguments which will override the defaults
#'   passed to [visreg::visreg()], The default arguments are \code{list(fit = x,
#'   trans = exp, plot = TRUE, rug = FALSE, alpha = 1, partial = FALSE, overlay
#'   = TRUE)}. For example, if you want a 95% confidence band, specify
#'   \code{hazard.params = list(alpha = 0.05)}. Note that The `cond` argument
#'   must be provided as a named list. Each element of that list specifies the
#'   value for one of the terms in the model; any elements left unspecified are
#'   filled in with the median/most common category. Only used for
#'   `type="hazard"`. All other argument are used for `type="hr"`. Note that the
#'   `visreg` package must be installed for `type="hazard"`.
#' @param newdata Required for \code{type="hr"}. The \code{newdata} argument is
#'   the "unexposed" group, while the exposed group is defined by either: (i) a
#'   change (defined by the \code{increment} argument) in a variable in newdata
#'   defined by the \code{var} argument ; or (ii) an exposed function that takes
#'   a data-frame and returns the "exposed" group (e.g. \code{exposed =
#'   function(data) transform(data, treat=1)}). This is a generalization of the
#'   behavior of the rstpm2 plot function. It allows both numeric and factor
#'   variables to be incremented or decremented. See references for rstpm2
#'   package. Only used for `type="hr"`
#' @param var specify the variable name for the exposed/unexposed (name is given
#'   as a character variable). If this argument is missing, then the
#'   \code{exposed} argument must be specified. This is the variable which will
#'   be incremented by the \code{increment} argument to give the exposed
#'   category. If \code{var} is coded as a factor variable, then
#'   \code{increment=1} will return the next level of the variable in
#'   \code{newdata}. \code{increment=2} will return two levels above, and so on.
#'   If the value supplied to \code{increment} is greater than the number of
#'   levels, this will simply return the max level. You can also decrement the
#'   categorical variable by specifying a negative value, e.g.,
#'   \code{increment=-1} will return one level lower than the value in
#'   \code{newdata}. If \code{var} is a numeric, than \code{increment} will
#'   increment (if positive) or decrement (if negative) by the supplied value.
#'   Only used for `type="hr"`.
#' @param increment Numeric value indicating how much to increment (if positive)
#'   or decrement (if negative) the \code{var} variable in \code{newdata}. See
#'   \code{var} argument for more details. Default is 1. Only used for
#'   `type="hr"`.
#' @param exposed function that takes \code{newdata} and returns the exposed
#'   dataset (e.g. function(data) transform(data, treat = 1)). This argument
#'   takes precedence over the \code{var} argument, i.e., if both \code{var} and
#'   \code{exposed} are correctly specified, only the \code{exposed} argument
#'   will be used. Only used for `type="hr"`.
#' @param xvar Variable to be used on x-axis for hazard ratio plots. If NULL,
#'   the function defaults to using the time variable used in the call to
#'   \code{fitSmoothHazard}. In general, this should be any continuous variable
#'   which has an interaction term with another variable. Only used for
#'   `type="hr"`.
#' @param ci Logical; if TRUE confidence bands are calculated. Only available
#'   for `family="glm"` and `family="gam"`, and only used for `type="hr"`,
#'   Default: !add. Confidence intervals for hazard ratios are calculated using
#'   the Delta Method.
#' @param ci.lvl Confidence level. Must be in (0,1), Default: 0.95. Only used
#'   for `type="hr"`.
#' @param ci.col Confidence band color. Only used if argument `ci=TRUE`,
#'   Default: 'grey'. Only used for `type="hr"`.
#' @param rug Logical. Adds a rug representation (1-d plot) of the event times
#'   (only for `status=1`), Default: !ci. Only used for `type="hr"`.
#' @return a plot of the hazard function or hazard ratio. For `type="hazard"`, a
#'   `data.frame` (returned invisibly) of the original data used in the fitting
#'   along with the data used to create the plots including `predictedhazard`
#'   which is the predicted hazard for a given covariate pattern and time.
#'   `predictedloghazard` is the predicted hazard on the log scale. `lowerbound`
#'   and `upperbound` are the lower and upper confidence interval bounds on the
#'   hazard scale (i.e. used to plot the confidence bands). `standarderror` is
#'   the standard error of the log hazard or log hazard ratio (only if
#'   `family="glm"` or `family="gam"`). For `type="hr"`, `log_hazard_ratio` and
#'   `hazard_ratio` is returned, and if `ci=TRUE`, `standarderror` (on the log
#'   scale) and `lowerbound` and `upperbound` of the `hazard_ratio` are
#'   returned.
#' @details This function has only been thoroughly tested for `family="glm"`. If
#'   the user wants more customized plot aesthetics, we recommend saving the
#'   results to a `data.frame` and using  the graphical package of their choice.
#' @examples
#' if (requireNamespace("splines", quietly = TRUE)) {
#' data("simdat") # from casebase package
#' library(splines)
#' simdat <- transform(simdat[sample(1:nrow(simdat), size = 200),],
#'                     treat = factor(trt, levels = 0:1,
#'                     labels = c("control","treatment")))
#'
#' fit_numeric_exposure <- fitSmoothHazard(status ~ trt*bs(eventtime),
#'                                         data = simdat,
#'                                         ratio = 1,
#'                                         time = "eventtime")
#'
#' fit_factor_exposure <- fitSmoothHazard(status ~ treat*bs(eventtime),
#'                                        data = simdat,
#'                                        ratio = 1,
#'                                        time = "eventtime")
#'
#' newtime <- quantile(fit_factor_exposure[["data"]][[fit_factor_exposure[["timeVar"]]]],
#'                     probs = seq(0.05, 0.95, 0.01))
#'
#' par(mfrow = c(1,3))
#' plot(fit_numeric_exposure,
#'      type = "hr",
#'      newdata = data.frame(trt = 0, eventtime = newtime),
#'      exposed = function(data) transform(data, trt = 1),
#'      xvar = "eventtime",
#'      ci = TRUE)
#'
#' #by default this will increment `var` by 1 for exposed category
#' plot(fit_factor_exposure,
#'      type = "hr",
#'      newdata = data.frame(treat = factor("control",
#'               levels = c("control","treatment")), eventtime = newtime),
#'      var = "treat",
#'      increment = 1,
#'      xvar = "eventtime",
#'      ci = TRUE,
#'      ci.col = "lightblue",
#'      xlab = "Time",
#'      main = "Hazard Ratio for Treatment",
#'      ylab = "Hazard Ratio",
#'      lty = 5,
#'      lwd = 7,
#'      rug = TRUE)
#'
#'
#' # we can also decrement `var` by 1 to give hazard ratio for control/treatment
#' result <- plot(fit_factor_exposure,
#'                type = "hr",
#'                newdata = data.frame(treat = factor("treatment",
#'                                     levels = c("control","treatment")),
#'                                     eventtime = newtime),
#'                var = "treat",
#'                increment = -1,
#'                xvar = "eventtime",
#'                ci = TRUE)
#'
#' # see data used to create plot
#' head(result)
#' }
#' @seealso \code{\link[utils]{modifyList}}, [casebase::fitSmoothHazard()],
#'   [graphics::par()], [visreg::visreg()]
#' @rdname plot.singleEventCB
#' @references Mark Clements and Xing-Rong Liu (2019). rstpm2: Smooth Survival
#'   Models, Including Generalized Survival Models. R package version 1.5.1.
#'   https://CRAN.R-project.org/package=rstpm2
#'
#'   Breheny P and Burchett W (2017). Visualization of Regression Models Using
#'   visreg. The R Journal, 9: 56-71.
#' @export
#' @importFrom utils modifyList
plot.singleEventCB <- function(x, ...,
                               type = c("hazard", "hr"),
                               hazard.params = list(),
                               newdata,
                               exposed,
                               increment = 1,
                               var,
                               xvar = NULL,
                               ci = FALSE,
                               ci.lvl = 0.95,
                               rug = !ci,
                               ci.col = "grey") {
    # Switch back call slot
    # otherwise we get an error from visreg
    x$call <- x$lower_call
    type <- match.arg(type)

    if (type == "hazard") {
        if (!requireNamespace("visreg", quietly = TRUE)) {
            stop("visreg package needed for this function. please install it first.")
        }

        tt <- do.call(visreg::visreg, utils::modifyList(
            list(
                fit = x,
                trans = exp,
                plot = TRUE,
                rug = FALSE,
                alpha = 1,
                partial = FALSE,
                overlay = TRUE,
                print.cond = TRUE
            ),
            hazard.params
        ))
    }


    if (type == "hr") {
        if (is.null(newdata) && type %in% c("hr")) {
            stop("Prediction using type 'hr' requires newdata to be specified.")
        }

        check_arguments_hazard(
            object = x, newdata = newdata, plot = FALSE,
            ci = ci, ci.lvl = ci.lvl
        )

        if (missing(exposed) & missing(var)) {
            stop("One of 'var' or 'exposed' arguments must be specified.")
        }

        if (missing(exposed)) {
            if (var %ni% names(newdata)) {
                stop(sprintf("%s not found in 'newdata'.", var))
            }
            newdata2 <- incrVar(var, increment = increment)(newdata)
        } else if (is.function(exposed)) {
            if (!missing(var)) warning("'var' argument is being ignored since 'exposure' argment has been correctly specified.")
            newdata2 <- exposed(newdata)
        } else if (!is.function(exposed) & !missing(var)) {
            if (var %ni% names(newdata)) {
                stop(sprintf("%s not found in 'newdata'.", var))
            }
            warning("'exposure' argument ignored since it is not correctly specified. Using 'var' argument instead.")
            newdata2 <- incrVar(var, increment = increment)(newdata)
        } else {
            stop("incorrect specification for 'exposed' argument. see help page for details.")
        }

        check_arguments_hazard(
            object = x, newdata = newdata2, plot = FALSE,
            ci = ci, ci.lvl = ci.lvl
        )

        tt <- plotHazardRatio(
            x = x, newdata = newdata, newdata2 = newdata2, ci = ci,
            ci.lvl = ci.lvl, ci.col = ci.col, rug = rug, xvar = xvar, ...
        )
    }
    invisible(tt)
}

#' @title Plot Cumulative Incidence and Survival Curves
#' @description Plot method for objects returned by the \code{absoluteRisk}
#'   function. Current plot types are cumulative incidence and survival
#'   functions.
#' @param x Fitted object of class `absRiskCB`. This is the result from the
#'   [casebase::absoluteRisk()] function.
#' @param ... further arguments passed to `matplot`. Only used if
#'   \code{gg=FALSE}.
#' @param xlab xaxis label, Default: 'time'
#' @param ylab yaxis label. By default, this will use the `"type"` attribute of
#'   the `absRiskCB` object
#' @param type Line type. Only used if `gg = FALSE`. This argument gets passed
#'   to [graphics::matplot()]. Default: 'l'
#' @param gg Logical for whether the `ggplot2` package should be used for
#'   plotting. Default: TRUE
#' @param id.names Optional character vector used as legend key when `gg=TRUE`.
#'   If missing, defaults to V1, V2, ...
#' @param legend.title Optional character vector of the legend title. Only used
#'   if `gg = FALSE`. Default is `'ID'`
#' @return A plot of the cumulative incidence or survival curve
#' @seealso \code{\link[graphics]{matplot}},
#'   \code{\link[casebase]{absoluteRisk}},
#'   \code{\link[data.table]{as.data.table}}, \code{\link[data.table]{setattr}},
#'   \code{\link[data.table]{melt.data.table}}
#' @rdname absoluteRisk
#' @export
#' @importFrom graphics matplot
#' @importFrom ggplot2 ggplot aes geom_line labs theme
#' @importFrom data.table as.data.table setnames melt
#' @examples
#' # Plot CI curves----
#' library(ggplot2)
#' data("brcancer")
#' mod_cb_tvc <- fitSmoothHazard(cens ~ estrec*log(time) +
#'                                 horTh +
#'                                 age +
#'                                 menostat +
#'                                 tsize +
#'                                 tgrade +
#'                                 pnodes +
#'                                 progrec,
#'                               data = brcancer,
#'                               time = "time", ratio = 1)
#' smooth_risk_brcancer <- absoluteRisk(object = mod_cb_tvc,
#'                                      newdata = brcancer[c(1,50),])
#'
#' class(smooth_risk_brcancer)
#' plot(smooth_risk_brcancer)
plot.absRiskCB <- function(x, ...,
                           xlab = "time",
                           ylab = ifelse(attr(x, "type") == "CI",
                                         "cumulative incidence",
                                         "survival probability"),
                           type = "l",
                           gg = TRUE,
                           id.names,
                           legend.title) {
    # output of absoluteRisk will always have a column named time
    # x equals linearRisk
    # ======================

    if (!gg) {
        graphics::matplot(x = x[, "time"],
                          y = x[, -which(colnames(x) == c("time"))],
                          type = type,
                          xlab = xlab,
                          ylab = ylab,
                          ...
        )

    } else {

        ID <- value <- NULL
        DT <- data.table::as.data.table(x)

        if (!missing(id.names)) {
            names_to_change <- grep("V", colnames(DT), value = TRUE)
            if (length(names_to_change) != length(id.names)) {
                warning("length of 'id.names' not equal to number of covariate profiles. ignoring 'id.names' argument")
            } else {
                data.table::setnames(DT, old = names_to_change, new = id.names)
            }
        }

        DT.m <- data.table::melt(DT, id.vars = "time", variable.name = "ID")

        ggplot(DT.m, aes(x = time, y = value, color = ID)) +
            geom_line() + theme_cb() +
            labs(color = ifelse(missing(legend.title), "ID", legend.title),
                 x = xlab, y = ylab) +
            theme(legend.position = "bottom")
    }

}
