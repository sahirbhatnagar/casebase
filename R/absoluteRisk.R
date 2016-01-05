# This is where all functions related to absolute risk estimation should appear

#' Compute absolute risks using the fitted hazard function.
#'
#' Using the output of the function \link{\code{fitSmoothHazard}}, we can
#' compute absolute risks by integrating the fitted hazard function over a time
#' period and then coverting this to an estimated survival for each individual.
#'
#' It is currently assumed that the time variable corresponds to the second, and
#' only the second, column of the design matrix. In particular, splines may not
#' be used at this time to fit a semi-parametric hazard function.
#'
#' @param object Output of function \link{\code{fitSmoothHazard}}.
#' @param time Upper bound of the interval over which to compute the absolute
#'   risk. It is assumed that \code{time} is given in the same units as the time
#'   variable in the dataset used to fit the hazard function.
#' @param newdata Optionally, a data frame in which to look for variables with
#'   which to predict. If omitted, the mean absolute risk is returned.
#' @param method Method used for integration. Defaults to \code{"quadrature"},
#'   which simply calls the function \link{\code{integrate}}. The only other
#'   option is \code{"montecarlo"}, which implements Monte-Carlo integration.
#' @param nsamp Maximal number of subdivisions (if \code{method = "quadrature"})
#'   or number of sampled points (if \code{method = "montecarlo}).
#' @return Returns the mean absolute risk at the user-supplied time.
absoluteRisk <- function(object, time, newdata = NULL, method = c("quadrature", "montecarlo"), nsamp=100) {
    method <- match.arg(method)
    # Create hazard function
    lambda <- function(x, pred, beta) {
        # Note: we don't include the offset in the estimation
        return(as.numeric(exp(beta[1] + beta[2] * x + crossprod(beta[3:length(beta)], pred))))
    }
    modelMat <- object$x[,-c(1, 2)]

    surv <- rep(NA, nobs)
    if (method == "quadrature") {
        for (i in 1:nobs) {
            surv[i] <- exp(-integrate(lambda, lower=0, upper=time, pred=modelMat[i,], beta=coefficients(object),
                                      subdivisions = nsamp)$value)
        }
    }
    if (method == "montecarlo") {
        sampledPoints <- runif(nsamp) * time
        for (i in 1:nobs) {
            surv[i] <- exp(-mean(lambda(sampledPoints, pred=modelMat[i,], beta=coefficients(object))))
        }
    }

    absRisk <- mean(1.0 - surv)
    return(absRisk)
}





