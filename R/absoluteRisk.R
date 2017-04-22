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
#' @param ... Extra parameters. Currently these are simply ignored.
#' @return Returns the estimated absolute risk for the user-supplied covariate profiles. This will
#'   be stored in a 2- or 3-dimensional array, depending on the input. See details.
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
#' out_linear <- fitSmoothHazard(event ~ time + z, DT)
#' out_log <- fitSmoothHazard(event ~ log(time) + z, DT)
#'
#' linear_risk <- absoluteRisk(out_linear, time = 10, newdata = data.table("z"=c(0,1)))
#' log_risk <- absoluteRisk(out_log, time = 10, newdata = data.table("z"=c(0,1)))
absoluteRisk <- function(object, ...) UseMethod("absoluteRisk")

#' @rdname absoluteRisk
#' @export
absoluteRisk.default <- function(object, ...) {
    stop("This function should be used with an object of class glm or CompRisk",
         call. = TRUE)
}

#' @rdname absoluteRisk
#' @export
absoluteRisk.glm <- function(object, time, newdata, method = c("montecarlo", "numerical"), nsamp = 1000, ...) {
    method <- match.arg(method)
    meanAR <- FALSE

    # Create hazard function
    lambda <- function(x, fit, newdata) {
        # Note: the offset should be set to zero when estimating the hazard.
        newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                               row.names = as.character(1:length(x)))
        newdata2[object$timeVar] <- x
        withCallingHandlers(pred <- predict(fit, newdata2),
                            warning = handler_bsplines)
        return(as.numeric(exp(pred)))
    }

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
        meanAR <- TRUE
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
                    output[i, j + 1] <- (time_ordered[i] - time_ordered[i - 1]) * mean(lambda(sampledPoints * (time_ordered[i] -
                                                                                                               time_ordered[i - 1]) +
                                                                                              time_ordered[i - 1], fit = object, newdata = newdata[j,,drop = FALSE]))
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

#' @rdname absoluteRisk
#' @export
absoluteRisk.CompRisk <- function(object, time, newdata, method = c("montecarlo", "numerical"), nsamp = 1000, ...) {
    # stop("absoluteRisk is not currently implemented for competing risks",
    #      call. = FALSE)
    method <- match.arg(method)
    meanAR <- FALSE

    if (missing(newdata)) {
        # Should we use the whole case-base dataset or the original one?
        if (is.null(object@originalData)) {
            stop("Can't estimate the mean absolute risk without the original data. See documentation.",
                 call. = FALSE)
        }
        newdata <- object@originalData
        # colnames(data)[colnames(data) == "event"] <- "status"
        # Next commented line will break on data.table
        # newdata <- newdata[, colnames(newdata) != "time"]
        unselectTime <- (names(newdata) != object@timeVar)
        newdata <- subset(newdata, select = unselectTime)
        meanAR <- TRUE
    }
    ###################################################
    # In competing risks, we can get a cumulative
    # incidence function using a nested double integral
    # f_j = lambda_j * Survival
    # F_j = P(T <= t, J = j : covariates) = int_0^t f_j
    ###################################################
    J <- length(object@typeEvents) - 1
    time_ordered <- unique(c(0, sort(time)))
    output <- array(NA, dim = c(length(time_ordered), nrow(newdata) + 1, J))
    output[,1,] <- time_ordered
    dimnames(output) <- list(rep("", length(time_ordered)),
                             c("time", rep("", nrow(newdata))),
                             paste("event", object@typeEvents[-1], sep = "="))
    output[1,-1,] <- 0

    # 1. Compute overall survival
    overallLambda <- function(x, object, newdata) {
        # Note: the offset should be set to zero when estimating the hazard.
        newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                               row.names = as.character(1:length(x)))
        newdata2[object@timeVar] <- x
        # predictvglm doesn't like offset = 0
        withCallingHandlers(pred <- VGAM::predictvglm(object, newdata2),
                            warning = handler_offset)
        return(as.numeric(rowSums(exp(pred))))
    }

    if (length(time_ordered) > 1) {

        if (method == "numerical") {
            overallSurv <- function(x, object, newdata) {
                exp(-integrate(overallLambda, lower = 0, upper = x[1], object = object, newdata = newdata,
                               subdivisions = nsamp)$value)
            }
            # 2. Compute individual subdensities f_j
            subdensities <- vector("list", length = J)
            subdensity_template <- function(x, object, newdata, index) {
                newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                                       row.names = as.character(1:length(x)))
                newdata2[object@timeVar] <- x
                # predictvglm doesn't like offset = 0
                withCallingHandlers(lambdas <- VGAM::predictvglm(object, newdata2),
                                    warning = handler_offset)
                exp(lambdas[,index]) * overallSurv(x, object = object, newdata2)
            }
            for (j in 1:J) {
                subdensities[[j]] <- partialize(subdensity_template, index = j)
            }

            # 3. Compute cumulative incidence functions F_j
            for (i in 1:nrow(newdata)) {
                for (j in 1:J) {
                    for (k in 2:length(time_ordered)) {
                        output[k,i + 1,j] <- integrate(subdensities[[j]], lower = time_ordered[k - 1],
                                                       upper = time_ordered[k],
                                                       object = object, newdata = newdata[i,,drop = FALSE],
                                                       subdivisions = nsamp)$value
                    }
                    output[,i,j] <- cumsum(output[,i,j])
                }
            }
        }

        if (method == "montecarlo") {
            sampledPoints <- runif(nsamp)
            overallSurv <- function(x, object, newdata) {
                sampledPoints <- runif(nsamp) * x
                exp(-x * mean(overallLambda(sampledPoints, object = object, newdata = newdata)))
            }
            for (i in 1:nrow(newdata)) {
                # output[i, ] <- time * colMeans(subdensity_mat(sampledPoints, object=object, newdata=newdata[i,]))
                for (k in 2:length(time_ordered)) {
                    x_vect <- sampledPoints * (time_ordered[k] - time_ordered[k - 1]) + time_ordered[k - 1]
                    surv <- overallSurv(x_vect, object, newdata[i,,drop = FALSE])
                    newdata2 <- data.frame(newdata[i,,drop = FALSE], offset = rep_len(0, length(x_vect)),
                                           row.names = as.character(1:length(x_vect)))
                    newdata2[object@timeVar] <- x_vect
                    withCallingHandlers(pred <- exp(VGAM::predictvglm(object, newdata2)),
                                        warning = handler_offset)

                    output[k, i + 1, ] <- (time_ordered[k] - time_ordered[k - 1]) * colMeans(surv * pred)
                }
                # if k==1, there was only one time point and we don't need to sum the contributions
                if (k != 1) output[,i,] <- apply(output[,i,], 2, cumsum)
            }

        }

    }

    # Reformat output when only one time point
    if (length(time) == 1) {
        if (time == 0) {
            output <- output[1,-1,,drop = FALSE]
        } else {
            output <- output[2,-1,,drop = FALSE]
        }
        dimnames(output)[[1]] <- time
    }

    # If there is only one time point, we should drop a dimension and return a matrix
    # output <- drop(output)

    return(output)
}
