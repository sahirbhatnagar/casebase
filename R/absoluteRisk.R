# This is where all functions related to absolute risk estimation should appear

#' Compute absolute risks using the fitted hazard function.
#'
#' Using the output of the function \code{fitSmoothHazard}, we can compute
#' absolute risks by integrating the fitted hazard function over a time period
#' and then converting this to an estimated survival for each individual.
#'
#' If \code{newdata = "typical"}, we create a typical covariate profile for the
#' absolute risk computation. This means that we take the median for numerical
#' and date variables, and we take the most common level for factor variables.
#'
#' In general, the output will include a column corresponding to the provided
#' time points. Some modifications of the \code{time} vector are done:
#' \code{time=0} is added, the time points are ordered, and duplicates are
#' removed. All these modifications simplify the computations and give an output
#' that can easily be used to plot risk curves.
#'
#' If there is no competing risk, the output is a matrix where each column
#' corresponds to the several covariate profiles, and where each row corresponds
#' to a time point. If there are competing risks, the output will be a
#' 3-dimensional array, with the third dimension corresponding to the different
#' events.
#'
#' The numerical method should be good enough in most situation, but Monte Carlo
#' integration can give more accurate results when the estimated hazard function
#' is not smooth (e.g. when modeling with time-varying covariates).
#'
#' @param object Output of function \code{\link{fitSmoothHazard}}.
#' @param time A vector of time points at which we should compute the absolute
#'   risks.
#' @param newdata Optionally, a data frame in which to look for variables with
#'   which to predict. If omitted, the mean absolute risk is returned.
#'   Alternatively, if \code{newdata = "typical"}, the absolute risk will be
#'   computed at a "typical" covariate profile (see Details).
#' @param method Method used for integration. Defaults to \code{"numerical"},
#'   which uses the trapezoidal rule to integrate over all time points together.
#'   The only other option is \code{"montecarlo"}, which implements Monte-Carlo
#'   integration.
#' @param nsamp Maximal number of subdivisions (if \code{method = "numerical"})
#'   or number of sampled points (if \code{method = "montecarlo"}).
#' @param n.trees Number of trees used in the prediction (for class \code{gbm}).
#' @param s Value of the penalty parameter lambda at which predictions are
#'   required (for class \code{cv.glmnet}).
#' @param onlyMain Logical. For competing risks, should we return absolute risks
#'   only for the main event of interest? Defaults to \code{TRUE}.
#' @param type Type of output. Can be \code{"CI"} (output is on the cumulative
#'   incidence scale) or \code{"survival"} (output is on the survival scale,
#'   i.e. 1- CI)
#' @param addZero Logical. Should we add time = 0 at the beginning of the
#'   output? Defaults to \code{TRUE}.
#' @param ntimes Number of time points (only used if \code{time} is missing).
#' @param ... Extra parameters. Currently these are simply ignored.
#' @return If \code{time} was provided, returns the estimated absolute risk for
#'   the user-supplied covariate profiles. This will be stored in a matrix or a
#'   higher dimensional array, depending on the input (see details). If both
#'   \code{time} and \code{newdata} are missing, returns the original data
#'   with a new column containing the risk estimate at failure times.
#' @export
#' @importFrom stats formula delete.response setNames
#' @examples
#' # Simulate censored survival data for two outcome types
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
#' DT <- data.table(z = rbinom(nobs, 1, 0.5))
#' DT[,`:=` ("t_event" = rweibull(nobs, 1, b1),
#'           "t_comp" = rweibull(nobs, 1, b2))]
#' DT[,`:=`("event" = 1 * (t_event < t_comp) + 2 * (t_event >= t_comp),
#'          "time" = pmin(t_event, t_comp))]
#' DT[time >= tlim, `:=`("event" = 0, "time" = tlim)]
#'
#' out_linear <- fitSmoothHazard(event ~ time + z, DT, ratio = 10)
#'
#' linear_risk <- absoluteRisk(out_linear, time = 10,
#'                             newdata = data.table("z" = c(0,1)))
absoluteRisk <- function(object, time, newdata,
                         method = c("numerical", "montecarlo"),
                         nsamp = 100, s = c("lambda.1se", "lambda.min"),
                         n.trees, onlyMain = TRUE,
                         type = c("CI", "survival"),
                         addZero = TRUE, ntimes = 100, ...) {
    if (!inherits(object, c("glm", "cv.glmnet", "gbm", "CompRisk"))) {
        stop(paste("object is of class", class(object)[1],
                   "\nabsoluteRisk should be used with an object of class glm,",
                   "cv.glmnet, gbm, or CompRisk"),
             call. = TRUE)
    }
    if (inherits(object, "CompRisk")) {
        return(absoluteRisk.CompRisk(object, time, newdata, method,
                                     nsamp = nsamp, onlyMain = onlyMain,
                                     type = type, addZero = addZero))
    }

    # Parse arguments
    method <- match.arg(method)
    type <- match.arg(type)
    if (is.numeric(s)) {
        s <- s[1]
        if (length(s) > 1) {
            warning(paste("More than one value for s has been",
                          "supplied. Only first entry will be used"))
            }
    } else if (is.character(s)) {
        s <- match.arg(s)
    }
    if (inherits(object, "gbm")) {
        if (missing(n.trees)) stop("n.trees is missing")
    } else n.trees <- NULL
    if (!missing(newdata) && is.character(newdata) && newdata == "typical") {
        newdata <- if (is.null(object$matrix.fit)) {
            get_typical(object$originalData)
        } else apply(object$originalData$x, 2, median, na.rm = TRUE)
    }

    # Call the correct function with correct parameters
    if (missing(newdata)) {
        if (missing(time)) {
            # If newdata and time are missing, compute risk for each subject
            # at their failure/censoring times
            return(estimate_risk(object, method, nsamp, s, n.trees, type = type,
                                 ...))
        } else {
            return(estimate_risk_newtime(object, time,
                                         method = method, nsamp = nsamp,
                                         s = s, n.trees = n.trees, type = type,
                                         addZero = addZero, ...))
        }
    } else {
        if (missing(time)) {
            # If only time is missing, compute risk for each row of newdata
            # at equidistant points
            max_time <- if (is.null(object$matrix.fit)) {
                max(object$originalData[object$timeVar][[1]])
            } else max(object$originalData$y[, object$timeVar])
            time <- seq(0, max_time, length.out = ntimes)
        }
        return(estimate_risk_newtime(object, time, newdata,
                                     method, nsamp, s, n.trees, type = type,
                                     addZero = addZero, ...))
    }
}

# Helper functions----
# estimate_risk computes one survival probability per covariate profile
# Currently, this is only called when we want the survival probabilities
# at failure times for the original data
estimate_risk <- function(object, method, nsamp, s, n.trees, type, ...) {
    newdata <- object$originalData
    if (inherits(newdata, "data.fit")) newdata <- newdata$x
    # Create risk variable and make sure it doesn't already exist
    riskVar <- ifelse(type == "CI", "risk", "survival")
    while (riskVar %in% names(newdata)) riskVar <- paste0(".", riskVar)
    # If both newdata and time are missing
    # compute risk at failure times
    time_vector <- if (is.null(object$matrix.fit)) {
        newdata[object$timeVar][[1]]
        } else object$originalData$y[, object$timeVar]
    # Create a vectorised hazard function
    if (inherits(object, "cv.glmnet") && !is.null(object$matrix.fit)) {
        hazard <- function(x, fit, newdata, s, n.trees, ...) {
            # Note: the offset should be set to zero when estimating the hazard.
            newdata_matrix <- newdata[, colnames(newdata) != fit$timeVar,
                                      drop = FALSE]
            # newdata_matrix matches output from fitSmoothHazard.fit
            temp_matrix <- model.matrix(update(fit$formula_time, ~ . - 1),
                                        setNames(data.frame(x), fit$timeVar))
            newdata_matrix <- as.matrix(cbind(temp_matrix,
                                              as.data.frame(newdata_matrix)))

            pred <- estimate_hazard(fit, newdata_matrix, plot = FALSE,
                                    ci = FALSE, s = s, n.trees = n.trees, ...)
            return(as.numeric(exp(pred)))
        }
    } else {
        hazard <- function(x, fit, newdata, s, n.trees, ...) {
            # Note: the offset should be set to zero when estimating the hazard.
            newdata2 <- data.frame(newdata)
            newdata2 <- newdata2[rep(1, length(x)), ]
            newdata2[fit$timeVar] <- x
            pred <- estimate_hazard(fit, newdata2, plot = FALSE, ci = FALSE,
                                    s = s, n.trees = n.trees, ...)
            return(as.numeric(exp(pred)))
        }
    }
    # Compute cumulative hazard at each failure time
    risk_res <- sapply(seq_len(nrow(newdata)), function(j) {
        integrate(hazard, lower = 0, upper = time_vector[j],
                  fit = object, subdivisions = nsamp,
                  newdata = newdata[j, , drop = FALSE],
                  s = s, n.trees = n.trees)$value
    })

    # Format the output depending on type
    if (is.data.frame(newdata)) {
        newdata[, riskVar] <- ifelse(type == "CI",
                                    1 - exp(-risk_res),
                                    exp(-risk_res))
    } else {
        newdata <- if (type == "CI") {
            cbind(1 - exp(-risk_res), newdata)
        } else cbind(exp(-risk_res), newdata)
        colnames(newdata)[1] <- riskVar
    }
    attr(newdata, "type") <- type
    return(newdata)
}

# estimate_risk_newtime computes the whole survival curve for each covariate
# profile
#' @importFrom data.table data.table :=
estimate_risk_newtime <- function(object, time, newdata, method, nsamp,
                                  s, n.trees, type, addZero, ...) {
    if (missing(newdata)) {
        # Should we use the whole case-base dataset or the original one?
        if (is.null(object$originalData)) {
            stop(paste("Cannot estimate the mean absolute risk without",
                       "the original data. See documentation."),
                 call. = FALSE)
        }
        newdata <- object$originalData
        unselectTime <- (names(newdata) != object$timeVar)
        newdata <- subset(newdata, select = unselectTime)
    }
    time_ordered <- unique(c(0, sort(time)))
    output <- matrix(NA, ncol = nrow(newdata) + 1, nrow = length(time_ordered))
    output[, 1] <- time_ordered
    colnames(output) <- c("time", rep("", nrow(newdata)))
    rownames(output) <- rep("", length(time_ordered))
    output[1, -1] <- 0

    if (length(time_ordered) > 1) {
        if (method == "numerical") {
            # Compute points at which we evaluate integral
            # Note: there's probably a more efficient way of choosing the knots
            knots <- seq(0, max(time_ordered),
                         length.out = (length(time_ordered) - 1) * nsamp)
            knots <- unique(sort(c(knots, time_ordered)))
            for (j in seq_len(nrow(newdata))) {
                # Extract current obs
                current_obs <- newdata[j, , drop = FALSE]
                # Use trapezoidal rule for integration----
                if (inherits(newdata, "matrix")) {
                    newdata2 <- current_obs[,
                                            colnames(newdata) != object$timeVar,
                                            drop = FALSE]
                    # newdata2 matches the output from fitSmoothHazard.fit
                    temp_matrix <- model.matrix(update(object$formula_time,
                                                       ~ . - 1),
                                                setNames(data.frame(knots),
                                                         object$timeVar))
                    newdata2 <- as.matrix(cbind(temp_matrix,
                                                as.data.frame(newdata2)))
                } else {
                    # Create data.table for prediction
                    newdata2 <- data.table::data.table(current_obs)
                    newdata2 <- newdata2[rep(1, length(knots))]
                    newdata2[, object$timeVar := knots]
                }
                pred <- estimate_hazard(object = object, newdata = newdata2,
                                        plot = FALSE, ci = FALSE,
                                        s = s, n.trees = n.trees, ...)
                # Compute integral using trapezoidal rule
                # First remove infinite values (e.g. with log(t))
                pred_exp <- exp(pred)
                pred_exp[which(pred_exp %in% c(Inf, -Inf))] <- 0
                which_knots <- knots %in% c(0, time_ordered)
                output[, j + 1] <- trap_int(knots, pred_exp)[which_knots]
            }

        }
        if (method == "montecarlo") {
            # Sample points at which we evaluate function
            knots <- runif(n = length(time_ordered) * nsamp,
                           min = 0, max = max(time_ordered))
            for (j in seq_len(nrow(newdata))) {
                # Extract current obs
                current_obs <- newdata[j, , drop = FALSE]
                # Create data.table for prediction
                newdata2 <- data.table(current_obs)
                newdata2 <- newdata2[rep(1, length(knots))]
                newdata2[, object$timeVar := knots]
                pred <- estimate_hazard(object = object, newdata = newdata2,
                                        plot = FALSE, ci = FALSE,
                                        s = s, n.trees = n.trees, ...)
                # Compute integral using MC integration
                pred_exp <- exp(pred)
                pred_exp[which(pred_exp %in% c(Inf, -Inf))] <- NA
                mean_values <- sapply(split(pred_exp,
                                            cut(knots, breaks = time_ordered)),
                                      mean, na.rm = TRUE)
                integral_estimates <- mean_values * diff(time_ordered)
                output[, j + 1] <- cumsum(c(0, integral_estimates))
            }
        }
        output[, -1] <- exp(-output[, -1])
        if (type == "CI") output[, -1] <- 1 - output[, -1]
    }
    # Sometimes montecarlo integration gives nonsensical probability estimates
    if (method == "montecarlo" && (any(output[, -1] < 0) |
                                   any(output[, -1] > 1))) {
        warning(paste("Some probabilities are out of range. Consider",
                      "increasing nsamp or using numerical integration"),
                call. = FALSE)
    }
    # Add time = 0?
    if (!addZero) output <- output[-1, , drop = FALSE]
    # Add class
    class(output) <- c("absRiskCB", class(output))
    attr(output, "type") <- type
    return(output)
}
