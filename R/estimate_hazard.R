estimate_hazard <- function(object, ...) UseMethod("estimate_hazard")

estimate_hazard.default <- function(object, ...) {
    stop(paste("This function should be used with an object of class",
               "glm, cv.glmnet, gbm"),
         call. = TRUE)
}

estimate_hazard.glm <- function(object, newdata, ci = FALSE, plot = FALSE,
                                ci.lvl = 0.95, ...) {
    check_arguments_hazard(object, newdata, plot, ci, ci.lvl)
    # Set offset to zero
    newdata$offset <- 0
    # Silence warnings about splines at t=0
    withCallingHandlers(pred <- predict(object, newdata, se.fit = ci,
                                        type = "link"),
                        warning = handler_bsplines)
    return(pred)
}

# The next function below makes the following assumption:
# Either time is already part of newdata (if it's a matrix)
# Or it's part of object$formula and will be added to newdata when we expand.
# Therefore, we need to make sure that assumption is verified before calling it.
estimate_hazard.cv.glmnet <- function(object, newdata, ci = FALSE, plot = FALSE,
                                      ci.lvl = 0.95, s = c("lambda.1se",
                                                           "lambda.min"), ...) {
    check_arguments_hazard(object, newdata, plot, ci, ci.lvl)
    if (is.numeric(s)) {
        if (length(s) > 1) {
            warning(paste("More than one value for s has been supplied.",
                          "Only first entry will be used"))
            }
        s <- s[1]
    } else if (is.character(s)) {
        s <- match.arg(s)
    }
    if (!inherits(newdata, "matrix")) {
        # Remove response variable because it won't be in newdata
        formula_pred <- formula(delete.response(terms(object$formula)))
        newdata_matrix <- prepareX(formula_pred, newdata)
    } else {
        newdata_matrix <- newdata
    }
    pred <- predict(object, newdata_matrix, s, newoffset = 0)

    return(pred)
}

estimate_hazard.gbm <- function(object, newdata, ci = FALSE, plot = FALSE,
                                ci.lvl = 0.95, n.trees, s = NULL, ...) {
    check_arguments_hazard(object, newdata, plot, ci, ci.lvl)
    # If gbm was fitted with an offset, predict.gbm ignores it
    # but still gives a warning. The following line silences this warning
    attr(object$Terms, "offset") <- NULL
    # Set offset to zero
    newdata$offset <- 0
    withCallingHandlers(pred <- predict(object, newdata, type = "link",
                                        n.trees, ...),
                        warning = handler_offset)
    return(pred)
}

check_arguments_hazard <- function(object, newdata, plot, ci, ci.lvl) {
    # For hazardPlot, we only want one row in newdata
    # for hazard ratio plot, we can have more than one row in newdata so
    # set plot = FALSE for hazard ratio plots in plot.singleEventCB function
    if (nrow(newdata) > 1 && plot) {
        newdata <- newdata[1, , drop = FALSE]
        warning(paste("More than 1 row supplied to 'newdata'.",
                      "Only the first row will be used."))
    }
    if (ci) {
        if (!data.table::between(ci.lvl, 0, 1, incbounds = FALSE))
            stop("ci.lvl must be between 0 and 1")
        if (!inherits(object, "glm")) {
            warning(sprintf(paste("Confidence intervals cannot be calculated",
                                  "for objects of class %s."),
                            class(object)[1]))
            ci <- FALSE
        }
        if (any(names(newdata) %in% c("standarderror", "lowerbound",
                                      "upperbound", "hazard_ratio",
                                      "log_hazard_ratio"))) {
            stop(paste("'standarderror','lowerbound', 'upperbound',",
                       "'hazard_ratio' and 'log_hazard_ratio' cannot be used",
                       "as column names in newdata. Rename it."))
        }
    } else {
        if (any(names(newdata) %in% c("hazard_ratio", "log_hazard_ratio")))
            stop(paste("'hazard_ratio' and 'log_hazard_ratio' cannot be used",
                       "as column names in newdata. Rename it."))
    }


    invisible(NULL)
}
