# This is where all functions related to absolute risk estimation should appear

#' Compute absolute risks using the fitted hazard function.
#'
#' Using the output of the function \code{fitSmoothHazard}, we can compute
#' absolute risks by integrating the fitted hazard function over a time period
#' and then coverting this to an estimated survival for each individual.
#'
#' In order to compute the mean absolute risk, the function \code{absoluteRisk}
#' needs the original dataset, i.e. the dataset before case-base sampling was
#' performed. This can be retrieved from the parameter \code{object} only if the
#' function \code{\link{fitSmoothHazard}} was run on the source dataset (and not
#' the output of \code{\link{sampleCaseBase}}). On the other hand, if the user
#' supplies the original dataset through the parameter \code{newdata}, the mean
#' absolute risk can be computed as the average of the output vector.
#'
#' The quadrature should be good enough in most situation, but Monte Carlo
#' integration can give more accurate results when the estimated hazard function
#' is not smooth (e.g. when modeling with time-varying covariates).
#'
#' @param object Output of function \code{\link{fitSmoothHazard}}.
#' @param time Upper bound of the interval over which to compute the absolute
#'   risk. It is assumed that \code{time} is given in the same units as the time
#'   variable in the dataset used to fit the hazard function.
#' @param newdata Optionally, a data frame in which to look for variables with
#'   which to predict. If omitted, the mean absolute risk is returned.
#' @param method Method used for integration. Defaults to \code{"quadrature"},
#'   which simply calls the function \code{\link{integrate}}. The only other
#'   option is \code{"montecarlo"}, which implements Monte-Carlo integration.
#' @param nsamp Maximal number of subdivisions (if \code{method = "quadrature"})
#'   or number of sampled points (if \code{method = "montecarlo"}).
#' @return Returns the mean absolute risk at the user-supplied time, if
#'   \code{newdata = NULL}, or the estimated absolute risk for the user-supplied
#'   covariate profiles.
#' @export
absoluteRisk <- function(object, ...) UseMethod("absoluteRisk")

#' @rdname absoluteRisk
#' @export
absoluteRisk.default <- function (object, ...) {
    stop("This function should be used with an object of class glm of compRisk",
         call. = FALSE)
}

#' @rdname absoluteRisk
#' @export
absoluteRisk.glm <- function(object, time, newdata = NULL, method = c("quadrature", "montecarlo"), nsamp=100) {
    method <- match.arg(method)
    meanAR <- FALSE
    # Create hazard function
    lambda <- function(x, fit, newdata) {
        # Note: the offset should be set to zero when estimating the hazard.
        newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                               row.names = as.character(1:length(x)))
        newdata2[object$timeVar] <- x
        return(as.numeric(exp(predict(fit, newdata2))))
    }

    if (is.null(newdata)) {
        # Should we use the whole case-base dataset or the original one?
        if(is.null(object$originalData)) {
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

    surv <- rep(NA, nrow(newdata))
    if (method == "quadrature") {
        for (i in 1:nrow(newdata)) {
            surv[i] <- exp(-integrate(lambda, lower=0, upper=time, fit=object, newdata=newdata[i,],
                                      subdivisions = nsamp)$value)
        }
    }
    if (method == "montecarlo") {
        sampledPoints <- runif(nsamp) * time
        for (i in 1:nrow(newdata)) {
            surv[i] <- exp(-mean(lambda(sampledPoints, fit=object, newdata=newdata[i,])))
        }
    }

    if (meanAR) {
        absRisk <- mean(1.0 - surv)
        return(absRisk)
    } else {
        absRisk <- 1.0 - surv
        names(absRisk) <- rownames(newdata)
        return(absRisk)
    }
}

#' @rdname absoluteRisk
#' @export
absoluteRisk.compRisk <- function(object, time, newdata = NULL, method = c("quadrature", "montecarlo"), nsamp=100) {
    # stop("absoluteRisk is not currently implemented for competing risks",
    #      call. = FALSE)
    method <- match.arg(method)
    meanAR <- FALSE

    if (is.null(newdata)) {
        # Should we use the whole case-base dataset or the original one?
        if(is.null(object$originalData)) {
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
    ###################################################
    # In competing risks, we can get a cumulative
    # incidence function using a nested double integral
    # f_j = lambda_j * Survival
    # F_j = P(T <= t, J = j : covariates) = int_0^t f_j
    ###################################################
    J <- length(object$typeEvents) - 1
    cumInc <- matrix(NA, nrow = nrow(newdata), ncol = J)

    # 1. Compute overall survival
    overallLambda <- function(x, fit, newdata) {
        # Note: the offset should be set to zero when estimating the hazard.
        newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                               row.names = as.character(1:length(x)))
        newdata2[object$timeVar] <- x
        # predictvglm doesn't like offset = 0
        withCallingHandlers(pred <- predictvglm(fit, newdata2),
                            warning = handler_offset)
        return(as.numeric(exp(rowSums(pred))))
    }

    if (method == "quadrature") {
        overallSurv <- function(time, fit, newdata) {
            exp(-integrate(overallLambda, lower=0, upper=time, fit=fit, newdata=newdata,
                           subdivisions = nsamp)$value)
        }
        # 2. Compute individual subdensities f_j
        subdensities <- vector("list", length = J)
        subdensity_template <- function(x, object, newdata, index) {
            newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                                   row.names = as.character(1:length(x)))
            newdata2[object$timeVar] <- x
            # predictvglm doesn't like offset = 0
            withCallingHandlers(lambdas <- predictvglm(object$model, newdata2),
                                warning = handler_offset)
            exp(lambdas[,index]) * overallSurv(x, fit = object$model, newdata2)
        }
        for (j in 1:J) {
            subdensities[[j]] <- partialize(subdensity_template, index = j, .lazy = FALSE)
        }

        # 3. Compute cumulative incidence functions F_j
        for (i in 1:nrow(newdata)) {
            for (j in 1:J) {
                cumInc[i,j] <- integrate(subdensities[[j]], lower = 0, upper = time,
                                         object=object, newdata=newdata[i,],
                                         subdivisions = nsamp)$value
            }
        }
    }

    if (method == "montecarlo") {
        overallSurv <- function(time, fit, newdata) {
            sampledPoints <- runif(nsamp) * time
            exp(-mean(overallLambda(sampledPoints, fit=fit, newdata=newdata)))
        }
        subdensity_mat <- function(x, object, newdata) {
            newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                                   row.names = as.character(1:length(x)))
            newdata2[object$timeVar] <- x
            # predictvglm doesn't like offset = 0
            withCallingHandlers(lambdas <- predictvglm(object$model, newdata2),
                                warning = handler_offset)
            exp(lambdas) * overallSurv(x, fit = object$model, newdata2)
        }
        sampledPoints <- runif(nsamp) * time
        for (i in 1:nrow(newdata)) {
            cumInc[i, ] <- colMeans(subdensity_mat(sampledPoints, object=object, newdata=newdata[i,]))
        }

    }

    return(cumInc)
}
