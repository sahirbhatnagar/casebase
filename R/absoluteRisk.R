# This is where all functions related to absolute risk estimation should appear

#' Compute absolute risks using the fitted hazard function.
#'
#' Using the output of the function \code{fitSmoothHazard}, we can compute absolute risks by
#' integrating the fitted hazard function over a time period and then converting this to an
#' estimated survival for each individual.
#'
#' If the user supplies the original dataset through the parameter \code{newdata}, the mean absolute
#' risk can be computed as the average of the output vector.
#'
#' In general, if \code{time} is a vector of length greater than one, the output will include a
#' column corresponding to the provided time points. Some modifications of the \code{time} vector
#' are done: \code{time=0} is added, the time points are ordered, and duplicates are removed. All
#' these modifications simplify the computations and give an output that can easily be used to plot
#' risk curves.
#'
#' On the other hand, if \code{time} corresponds to a single time point, the output does not include
#' a column corresponding to time.
#'
#' If there is no competing risk, the output is a matrix where each column corresponds to the
#' several covariate profiles, and where each row corresponds to a time point. If there are
#' competing risks, the output will be a 3-dimensional array, with the third dimension corresponding
#' to the different events.
#'
#' The numerical method should be good enough in most situation, but Monte Carlo integration can
#' give more accurate results when the estimated hazard function is not smooth (e.g. when modeling
#' with time-varying covariates). However, if there are competing risks, we strongly encourage the
#' user to select Monte-Carlo integration, which is much faster than the numerical method. (This is
#' due to the current implementation of the numerical method, and it may be improved in future
#' versions.)
#'
#' @param object Output of function \code{\link{fitSmoothHazard}}.
#' @param time A vector of time points at which we should compute the absolute risks.
#' @param newdata Optionally, a data frame in which to look for variables with which to predict. If
#'   omitted, the mean absolute risk is returned.
#' @param method Method used for integration. Defaults to \code{"montecarlo"}, which implements
#'   Monte-Carlo integration. The only other option is \code{"numerical"}, which simply calls the
#'   function \code{\link{integrate}}.
#' @param nsamp Maximal number of subdivisions (if \code{method = "numerical"}) or number of sampled
#'   points (if \code{method = "montecarlo"}).
#' @param n.trees Number of trees used in the prediction (for class \code{gbm}).
#' @param s Value of the penalty parameter lambda at which predictions are required (for class
#'   \code{cv.glmnet}).
#' @param onlyMain Logical. For competing risks, should we return absolute risks only for the main
#'   event of interest? Defaults to \code{TRUE}.
#' @param ... Extra parameters. Currently these are simply ignored.
#' @return If \code{time} was provided, returns the estimated absolute risk for the user-supplied
#'   covariate profiles. This will be stored in a 2- or 3-dimensional array, depending on the input
#'   (see details). If both \code{time} and \code{newdata} were provided, returns the original data
#'   with a new column containing the risk estimate at failure time.
#' @export
#' @examples
#' # Simulate censored survival data for two outcome types from exponential distributions
#' library(data.table)
#' set.seed(12345)
#' nobs <- 1000
#' tlim <- 20
#'
#' # simulation parameters
#' b1 <- 200
#' b2 <- 50
#'
#' # event type 0-censored, 1-event of interest, 2-competing event
#' # t observed time/endpoint
#' # z is a binary covariate
#' DT <- data.table(z=rbinom(nobs, 1, 0.5))
#' DT[,`:=` ("t_event" = rweibull(nobs, 1, b1),
#'           "t_comp" = rweibull(nobs, 1, b2))]
#' DT[,`:=`("event" = 1 * (t_event < t_comp) + 2 * (t_event >= t_comp),
#'          "time" = pmin(t_event, t_comp))]
#' DT[time >= tlim, `:=`("event" = 0, "time" = tlim)]
#'
#' out_linear <- fitSmoothHazard(event ~ time + z, DT, ratio = 10)
#' out_log <- fitSmoothHazard(event ~ log(time) + z, DT, ratio = 10)
#'
#' linear_risk <- absoluteRisk(out_linear, time = 10, newdata = data.table("z"=c(0,1)))
#' log_risk <- absoluteRisk(out_log, time = 10, newdata = data.table("z"=c(0,1)))
absoluteRisk <- function(object, ...) UseMethod("absoluteRisk")

#' @rdname absoluteRisk
#' @export
absoluteRisk.default <- function(object, ...) {
    stop("This function should be used with an object of class glm, cv.glmnet, gbm, or CompRisk",
         call. = TRUE)
}

#' @rdname absoluteRisk
#' @export
absoluteRisk.glm <- function(object, time, newdata, method = c("montecarlo", "numerical"), nsamp = 1000, ...) {
    method <- match.arg(method)

    # Create hazard function
    lambda <- function(x, fit, newdata) {
        # Note: the offset should be set to zero when estimating the hazard.
        newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                               row.names = as.character(1:length(x)))
        newdata2[fit$timeVar] <- x
        withCallingHandlers(pred <- predict(fit, newdata2),
                            warning = handler_bsplines)
        return(as.numeric(exp(pred)))
    }

    return(call_correct_estimate_risk_fn(lambda, object, time, newdata, method, nsamp))
    # return(estimate_risk(lambda, object, time, newdata, method, nsamp))
}

#' @rdname absoluteRisk
#' @export
absoluteRisk.gbm <- function(object, time, newdata, method = c("montecarlo", "numerical"), nsamp = 1000, n.trees, ...) {
    method <- match.arg(method)

    # Create hazard function
    lambda <- function(x, fit, newdata) {
        # Note: the offset should be set to zero when estimating the hazard.
        newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                               row.names = as.character(1:length(x)))
        newdata2[fit$timeVar] <- x
        withCallingHandlers(pred <- predict(fit, newdata2, n.trees, ...),
                            warning = handler_offset)
        return(as.numeric(exp(pred)))
    }

    return(call_correct_estimate_risk_fn(lambda, object, time, newdata, method, nsamp))
    # return(estimate_risk(lambda, object, time, newdata, method, nsamp))
}

#' @rdname absoluteRisk
#' @importFrom stats formula delete.response setNames
#' @export
absoluteRisk.cv.glmnet <- function(object, time, newdata, method = c("montecarlo", "numerical"),
                                   nsamp = 1000, s = c("lambda.1se","lambda.min"), ...) {
    method <- match.arg(method)
    if (is.numeric(s))
        s <- s[1]
    else if (is.character(s)) {
        s <- match.arg(s)
    }

    # Create hazard function
    if (is.null(object$matrix.fit)) {
        lambda <- function(x, fit, newdata) {
            # Note: the offset should be set to zero when estimating the hazard.
            newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                                   row.names = as.character(1:length(x)))
            newdata2[fit$timeVar] <- x
            formula_pred <- update(formula(delete.response(terms(fit$formula))), ~ . -1)
            newdata_matrix <- model.matrix(formula_pred, newdata2)
            pred <- predict(fit, newdata_matrix, s, newoffset = 0)
            return(as.numeric(exp(pred)))
        }
    } else {
        lambda <- function(x, fit, newdata) {
            # Note: the offset should be set to zero when estimating the hazard.
            # newdata_matrix <- cbind(x, newdata)
            newdata_matrix <- newdata[,colnames(newdata) != fit$timeVar, drop = FALSE]
            newdata_matrix <- as.matrix(cbind(as.data.frame(newdata_matrix),
                                              model.matrix(update(fit$formula_time, ~ . -1),
                                                           setNames(data.frame(x), fit$timeVar))))
            pred <- predict(fit, newdata_matrix, s, newoffset = 0)
            return(as.numeric(exp(pred)))
        }
    }

    return(call_correct_estimate_risk_fn(lambda, object, time, newdata, method, nsamp))
}

# Helper functions----
call_correct_estimate_risk_fn <- function(lambda, object, time, newdata, method, nsamp) {
    # Call the correct function with correct parameters
    if (missing(newdata)) {
        if (missing(time)) {
            return(estimate_risk(lambda, object))
        } else {
            return(estimate_risk_newtime(lambda, object, time,
                                         method = method, nsamp = nsamp))
        }
    } else {
        return(estimate_risk_newtime(lambda, object, time,
                                     newdata, method, nsamp))
    }
}

# The absolute risk methods create the lambda function and pass it to
# estimate_risk or estimate_risk_newtime
estimate_risk <- function(lambda, object) {
    newdata <- object$originalData
    if (inherits(newdata, "data.fit")) newdata <- newdata$x
    # If both newdata and time are missing
    # compute risk at failure times
    riskVar <- "risk"
    while (riskVar %in% names(newdata)) riskVar <- paste0(".", riskVar)

    time_vector <- if (is.null(object$matrix.fit)) {
        newdata[object$timeVar][[1]]
        } else object$originalData$y[,object$timeVar]
    risk_res <- sapply(seq_len(nrow(newdata)), function(j) {
        integrate(lambda, lower = 0, upper = time_vector[j],
                  fit = object, newdata = newdata[j,,drop = FALSE])$value
    })
    if (is.data.frame(newdata)) {
        newdata[,riskVar] <- 1 - exp(-risk_res)
    } else {
        newdata <- cbind(1 - exp(-risk_res), newdata)
        colnames(newdata)[1] <- riskVar
    }
    return(newdata)
}

estimate_risk_newtime <- function(lambda, object, time, newdata, method, nsamp) {
    if (missing(newdata)) {
        # Should we use the whole case-base dataset or the original one?
        if (is.null(object$originalData)) {
            stop("Can't estimate the mean absolute risk without the original data. See documentation.",
                 call. = FALSE)
        }
        newdata <- object$originalData
        # colnames(data)[colnames(data) == "event"] <- "status"
        # Next commented line will break on data.table
        # newdata <- newdata[, colnames(newdata) != "time"]
        unselectTime <- (names(newdata) != object$timeVar)
        newdata <- subset(newdata, select = unselectTime)
    }
    time_ordered <- unique(c(0, sort(time)))
    output <- matrix(NA, ncol = nrow(newdata) + 1, nrow = length(time_ordered))
    output[,1] <- time_ordered
    colnames(output) <- c("time", rep("", nrow(newdata)))
    rownames(output) <- rep("", length(time_ordered))
    output[1,-1] <- 0

    if (length(time_ordered) > 1) {
        # output[1,-1] <- 0
        if (method == "numerical") {
            for (j in 1:nrow(newdata)) {
                for (i in 2:length(time_ordered)) {

                    output[i, j + 1] <- integrate(lambda, lower = time_ordered[i - 1], upper = time_ordered[i],
                                                  fit = object, newdata = newdata[j,,drop = FALSE],
                                                  subdivisions = nsamp)$value
                }
                output[,j + 1] <- cumsum(output[,j + 1])
            }

        }
        if (method == "montecarlo") {
            sampledPoints <- runif(nsamp)
            for (j in 1:nrow(newdata)) {
                for (i in 2:length(time_ordered)) {
                    # output[i, j + 1] <- (time_ordered[i] - time_ordered[i - 1]) * mean(lambda(sampledPoints * (time_ordered[i] -
                    #                                                                                                time_ordered[i - 1]) +
                    #                                                                               time_ordered[i - 1], fit = object, newdata = newdata[j,,drop = FALSE]))
                    output[i, j + 1] <- integrate_mc(lambda, lower = time_ordered[i - 1], upper = time_ordered[i],
                                                     fit = object, newdata = newdata[j,,drop = FALSE],
                                                     subdivisions = nsamp)
                }
                output[,j + 1] <- cumsum(output[,j + 1])
            }
        }
        output[,-1] <- exp(-output[,-1])
        output[,-1] <- 1 - output[,-1]
    }

    # Reformat output when only one time point
    if (length(time) == 1) {
        if (time == 0) {
            output <- output[1,-1,drop = FALSE]
        } else {
            output <- output[2,-1,drop = FALSE]
        }
        rownames(output) <- time
    }
    return(output)
}
