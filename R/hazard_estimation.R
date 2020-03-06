estimate_hazard <- function(object, newdata, ci = FALSE, plot = FALSE, ...) {
    # For hazardPlot, we only want one row in newdata
    if (nrow(newdata) > 1 && plot) {
        newdata <- newdata[1, , drop = FALSE]
        warning("More than 1 row supplied to 'newdata'. Only the first row will be used.")
    }
    UseMethod("estimate_hazard")
}

estimate_hazard.default <- function(object, newdata, ci, plot, ...) {
    stop("This function should be used with an object of class glm, cv.glmnet, gbm",
         call. = TRUE)
}

estimate_hazard.glm <- function(object, newdata, ci, plot, ...) {
    withCallingHandlers(pred <- predict(object, newdata, type = "link"),
                        warning = handler_bsplines)
    return(pred)
}

estimate_hazard.cv.glmnet <- function(object, newdata, ci, plot,
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

estimate_hazard.gbm <- function(object, newdata, ci, plot, n.trees, ...) {
    withCallingHandlers(pred <- predict(object, newdata, type = "link",
                                        n.trees, ...),
                        warning = handler_offset)
    return(pred)
}
