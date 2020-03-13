estimate_hazard <- function(object, newdata, ci = FALSE, plot = FALSE, ci.lvl = 0.95, ...) {
    # For hazardPlot, we only want one row in newdata
    if (nrow(newdata) > 1 && plot) {
        newdata <- newdata[1, , drop = FALSE]
        warning("More than 1 row supplied to 'newdata'. Only the first row will be used.")
    }
    if (ci) {
        if (!data.table::between(ci.lvl, 0,1, incbounds = FALSE))
            stop("ci.lvl must be between 0 and 1")
        if (!inherits(object, "glm")) {
            warning(sprintf("Confidence intervals cannot be calculated for objects of class %s.",
                            class(object)[1]))
            ci <- FALSE
        }
        if (any(names(newdata) %in% c("standarderror","lowerbound","upperbound")))
            stop("'standarderror','lowerbound' and 'upperbound' cannot be used as column names in newdata. rename it.")
    }
    UseMethod("estimate_hazard")
}

estimate_hazard.default <- function(object, newdata, ci = FALSE, plot = FALSE, ci.lvl = 0.95, ...) {
    stop("This function should be used with an object of class glm, cv.glmnet, gbm",
         call. = TRUE)
}

estimate_hazard.glm <- function(object, newdata, ci = FALSE, plot = FALSE, ci.lvl = 0.95, ...) {
    # Set offset to zero
    newdata$offset <- 0
    # Silence warnings about splines at t=0
    withCallingHandlers(pred <- predict(object, newdata, se.fit = ci, type = "link"),
                        warning = handler_bsplines)
    return(pred)
}

estimate_hazard.cv.glmnet <- function(object, newdata, ci = FALSE, plot = FALSE, ci.lvl = 0.95,
                                   s = c("lambda.1se","lambda.min"), ...) {
    if (is.numeric(s))
        s <- s[1]
    else if (is.character(s)) {
        s <- match.arg(s)
    }
    if (!inherits(newdata, "matrix")) {
        # Remove the intercept to get the right design matrix
        formula_pred <- update(formula(delete.response(terms(object$formula))), ~ . -1)
        newdata_matrix <- model.matrix(formula_pred, newdata)
    } else {
        newdata_matrix <- newdata
    }
    pred <- predict(object, newdata_matrix, s, newoffset = 0)

    return(pred)
}

estimate_hazard.gbm <- function(object, newdata, ci = FALSE, plot = FALSE, ci.lvl = 0.95, n.trees, ...) {
    # If gbm was fitted with an offset, predict.gbm ignores it but still gives a warning
    # The following line silences this warning
    attr(object$Terms, "offset") <- NULL
    # Set offset to zero
    newdata$offset <- 0
    withCallingHandlers(pred <- predict(object, newdata, type = "link",
                                        n.trees, ...),
                        warning = handler_offset)
    return(pred)
}
