#' @title Plot Fitted Hazard Curve as a Function of Time
#' @description Visualize estimated hazard curves as a function of time with
#'   confidence intervals. This function takes as input, the result from the
#'   [casebase::fitSmoothHazard()] function. The user can also specify a
#'   sequence of times at which to estimate the hazard function. These plots are
#'   useful to visualize the non-proportional hazards, i.e., time dependent
#'   interactions with a covariate.
#' @param object Fitted object of class `glm`, `gam`, `cv.glmnet` or `gbm`. This
#'   is the result from the [casebase::fitSmoothHazard()] function.
#' @param newdata A data frame in which to look for variables with which to
#'   predict. This is required and must contain all the variables used in the
#'   model. Only one covariate profile can be used. If more than one row is
#'   provided, only the first row will be used.
#' @param type Type of plot. Currently, only "hazard" has been implemented.
#'   Default: c("hazard")
#' @param xlab x-axis label. Default: the name of the time variable from the
#'   fitted `object`.
#' @param breaks Number of points at which to estimate the hazard. This argument
#'   is only used if argument `times=NULL`. This function will calculate a
#'   sequence of times between the minimum and maximum of observed event times.
#'   Default: 100.
#' @param ci.lvl Confidence level. Must be in (0,1), Default: 0.95
#' @param ylab y-axis label. Default: NULL which means the function will put
#'   sensible defaults.
#' @param line.col Line color, Default: 1. See [graphics::par()] for details.
#' @param ci.col Confidence band color. Only used if argument `ci=TRUE`,
#'   Default: 'grey'
#' @param lty Line type. See [graphics::par()] for details, Default: par("lty")
#' @param add Logical; if TRUE add to an already existing plot; Default: FALSE
#' @param ci Logical; if TRUE confidence bands are calculated. Only available
#'   for `family="glm"` and `family="gam"`, Default: !add
#' @param rug Logical. Adds a rug representation (1-d plot) of the event times
#'   (only for `status=1`), Default: !add
#' @param s Value of the penalty parameter lambda at which predictions are
#'   required (for class \code{cv.glmnet} only). Only the first entry will be
#'   used if more than one numeric value is provided, Default: c("lambda.1se",
#'   "lambda.min")
#' @param times Vector of numeric values at which the hazard should be
#'   calculated. Default: NULL which means this function will use the minimum
#'   and maximum of observed event times with the `breaks` argument.
#' @param ... further arguments passed to [graphics::matplot()]
#' @details This is an earlier version of a function to plot hazards. We
#'   recommend instead using the plot method for objects returned by
#'   [casebase::fitSmoothHazard()]. See [casebase::plot.singleEventCB()].
#' @return a plot of the hazard function and a data.frame of original data used
#'   in the fitting along with the data used to create the plots including
#'   `predictedhazard` which is the predicted hazard for a given covariate
#'   pattern and time `predictedloghazard` is the predicted hazard on the log
#'   scale. `lowerbound` and `upperbound` are the lower and upper confidence
#'   interval bounds on the hazard scale (i.e. used to plot the confidence
#'   bands). `standarderror` is the standard error of the log hazard (only if
#'   `family="glm"` or `family="gam"`)
#' @seealso [casebase::fitSmoothHazard()]
#' @examples
#' data("simdat")
#' mod_cb <- fitSmoothHazard(status ~ trt * eventtime,
#'                                     time = "eventtime",
#'                                     data = simdat[1:200,],
#'                                     ratio = 1,
#'                                     family = "glm")
#'
#' results0 <- hazardPlot(object = mod_cb, newdata = data.frame(trt = 0),
#'            ci.lvl = 0.95, ci = FALSE, lty = 1, line.col = 1, lwd = 2)
#' head(results0)
#' hazardPlot(object = mod_cb, newdata = data.frame(trt = 1), ci = FALSE,
#'            ci.lvl = 0.95, add = TRUE, lty = 2, line.col = 2, lwd = 2)
#' legend("topleft", c("trt=0","trt=1"),lty=1:2,col=1:2,bty="y", lwd = 2)
#' @rdname hazardPlot
#' @export
#' @importFrom graphics lines matplot par polygon
#' @importFrom data.table between
#' @importFrom stats qnorm
hazardPlot <- function(object, newdata, type = c("hazard"), xlab = NULL,
                       breaks = 100, ci.lvl = 0.95, ylab = NULL, line.col = 1,
                       ci.col = "grey", lty = par("lty"), add = FALSE,
                       ci = !add, rug = !add, s = c("lambda.1se", "lambda.min"),
                       times = NULL, ...) {

    type <- match.arg(type)

    if (is.null(newdata))
        stop("newdata argument needs to be specified")

    if (!is.data.frame(newdata))
        stop("newdata must be a data.frame")

    if (nrow(newdata) > 1) {
        newdata <- newdata[1, , drop = FALSE]
        warning("More than 1 row supplied to 'newdata'. Only the first row will be used.")
    }

    if (any(names(newdata) %in% c("sequence_of_times", "predictedloghazard")))
        stop("sequence_of_times or predictedloghazard cannot be used as a colunm name in newdata. rename them to something else.")

    if (any(names(newdata) %in% object[["timeVar"]]))
        stop("%s cannot be used as a colunm name in newdata. remove it.",
             object[["timeVar"]])

    obj_class <- class(object)[2]

    if (!inherits(object, c("glm", "gam", "gbm", "cv.glmnet")))
        stop("object must be of class glm, gam, gbm or cv.glmnet")

    if (ci) {
        if (!data.table::between(ci.lvl, 0, 1, incbounds = FALSE))
            stop("ci.lvl must be between 0 and 1")
        if (!inherits(object, c("glm", "gam"))) {
            warning(sprintf("Confidence intervals cannot be calculated for objects of class %s.",
                            obj_class))
            ci <- FALSE
        }
        if (any(names(newdata) %in% c("standarderror", "lowerbound",
                                      "upperbound")))
            stop("'standarderror','lowerbound' and 'upperbound' cannot be used as column names in newdata. rename it.")
    }

    if (is.null(times)) {
        times <- object[["originalData"]][[object[["timeVar"]]]]
        times <- seq(min(times), max(times), length.out = breaks)
    } else {
        times <- times[order(times)]
    }

    newdata <- newdata[rep(seq_len(nrow(newdata)), each = length(times)), ,
                       drop = FALSE]
    newdata$sequence_of_times <- times
    names(newdata)[names(newdata) == "sequence_of_times"] <- object[["timeVar"]]

    newdata$offset <- 0
    # If gbm was fitted with an offset, predict.gbm ignores it
    # but still gives a warning. The following line silences this warning
    if (!inherits(object, "gbm")) attr(object$Terms, "offset") <- NULL

    preds <- switch(obj_class,
                    glm = {
                        pp <- predict(object, newdata = newdata,
                                      se.fit = ci, type = "link")
                        newdata$predictedloghazard <- if (ci) pp[["fit"]] else pp
                        pp
                    },
                    gam = {
                        pp <- predict(object, newdata = newdata,
                                      se.fit = ci, type = "link")
                        newdata$predictedloghazard <- if (ci) pp[["fit"]] else pp
                        pp
                     },
                    cv.glmnet = {
                        if (is.numeric(s)) {
                            if (length(s) > 1)  {
                                warning(paste("More than one value for s has",
                                              "been supplied. Only first entry",
                                              "will be used"))
                            }
                            s <- s[1]
                        } else if (is.character(s)) {
                            s <- match.arg(s)
                        }
                        # First extract RHS, then remove intercept
                        newx <- model.matrix(update(object$formula[-2],
                                                    ~ . - 1),
                                             newdata)
                        # the newoffset = 0 isn't required for now as glmnet
                        # is not using the offset because of the hack
                        pp <- predict(object, newx = newx, s = s, newoffset = 0)
                        newdata$predictedloghazard <- pp
                        pp
                    },
                    gbm = {
                        # If gbm was fitted with an offset
                        # predict.gbm ignores it but still gives a warning
                        # The following line silences this warning
                        attr(object$Terms, "offset") <- NULL
                        pp <- predict(object, newdata, n.trees = object$n.trees)
                        newdata$predictedloghazard <- pp
                        pp
                    }
    )

    if (ci) {
        std_err <- newdata$standarderror <- preds$se.fit
        newdata$lowerbound <- exp(newdata$predictedloghazard +
                                      qnorm(0.5 * (1 - ci.lvl)) * std_err)
        newdata$upperbound <- exp(newdata$predictedloghazard +
                                      qnorm(1 - 0.5 * (1 - ci.lvl)) * std_err)
    }

    newdata$predictedhazard <- exp(newdata$predictedloghazard)

    if (is.null(xlab))
        xlab <- object[["timeVar"]]
    if (is.null(ylab))
        ylab <- switch(type, hr = "Hazard ratio", hazard = "Hazard",
                       surv = "Survival", density = "Density",
                       sdiff = "Survival difference",
                       hdiff = "Hazard difference", marghaz = "Marginal hazard",
                       loghazard = "log(hazard)", link = "Linear predictor",
                       meansurv = "Mean survival", cumhaz = "Cumulative hazard",
                       meansurvdiff = "Difference in mean survival",
                       meanhr = "Mean hazard ratio", odds = "Odds",
                       or = "Odds ratio", margsurv = "Marginal survival",
                       marghr = "Marginal hazard ratio", haz = "Hazard",
                       fail = "Failure", meanhaz = "Mean hazard",
                       margfail = "Marginal failure",
                       af = "Attributable fraction",
                       meanmargsurv = "Mean marginal survival",
                       uncured = "Uncured distribution")

    ylims <- if (ci) {
        range(newdata$lowerbound, newdata$upperbound)
        } else range(newdata$predictedhazard)

    if (!add)
        matplot(newdata[[object[["timeVar"]]]], newdata[["predictedhazard"]],
                type = "n", xlab = xlab, ylab = ylab, ylim = ylims, ...)
    if (ci) {
        polygon(c(newdata[[object[["timeVar"]]]],
                  rev(newdata[[object[["timeVar"]]]])),
                c(newdata[["lowerbound"]], rev(newdata[["upperbound"]])),
                col = ci.col,
                border = ci.col)
        lines(newdata[[object[["timeVar"]]]], newdata[["predictedhazard"]],
              col = line.col, lty = lty, ...)
    } else lines(newdata[[object[["timeVar"]]]], newdata[["predictedhazard"]],
                 col = line.col, lty = lty, ...)
    if (rug) {
        events <- object[["originalData"]][[object[["eventVar"]]]]
        rug(object[["originalData"]][which(events == 1), ,
                                     drop = FALSE][[object[["timeVar"]]]],
            col = line.col, quiet = TRUE) # Silence warning about clipped values
    }
    return(invisible(newdata))
}
