# This is where all functions related to absolute risk estimation should appear

#' Compute absolute risks using the fitted hazard function.
#'
#' Using the output of the function \code{fitSmoothHazard}, we can
#' compute absolute risks by integrating the fitted hazard function over a time
#' period and then coverting this to an estimated survival for each individual.
#'
#' In order to compute the mean absolute risk, the function \code{absoluteRisk}
#' needs the original dataset, i.e. the dataset before case-base sampling was
#' performed. This can be retrieved from the parameter \code{object} only if the
#' function \code{\link{fitSmoothHazard}} was run on the source dataset (and not
#' the output of \code{\link{sampleCaseBase}}). On the other hand, if the user
#' supplies the original dataset through the parameter \code{newdata}, the mean
#' absolute risk can be computed as the average of the output vector.
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
absoluteRisk <- function(object, time, newdata = NULL, method = c("quadrature", "montecarlo"), nsamp=100) {
    method <- match.arg(method)
    meanAR <- FALSE
    # Create hazard function
    lambda <- function(x, fit, newdata) {
        # Note: the offset should be set to zero when estimating the hazard.
        newdata2 <- data.frame(time = x, newdata, offset = rep_len(0, length(x)),
                               row.names = as.character(1:length(x)))
        return(as.numeric(exp(predict(fit, newdata2))))
    }

    if (is.null(newdata)) {
        # Should we use the whole case-base dataset or the original one?
        if(is.null(object$originalData)) {
            stop("Can't estimate the mean absolute risk without the original data. See documentation.")
        }
        newdata <- object$originalData
        # colnames(data)[colnames(data) == "event"] <- "status"
        # Next commented line will break on data.table
        # newdata <- newdata[, colnames(newdata) != "time"]
        newdata <- subset(newdata, select = (colnames(newdata) != "time"))
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





