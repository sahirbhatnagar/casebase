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
#' @return Returns the mean absolute risk at the user-supplied time, if
#'   \code{newdata = NULL}, or the estimated absolute risk for the user-supplied
#'   covariate profiles.
absoluteRisk <- function(object, time, newdata = NULL, method = c("quadrature", "montecarlo"), nsamp=100) {
    method <- match.arg(method)
    # Create hazard function
    lambda <- function(x, fit, newdata) {
        # Note: the offset should be set to zero when estimating the hazard.
        newdata2 <- data.frame(time = x, newdata, offset = rep_len(0, length(x)),
                               row.names = as.character(1:length(x)))
        return(as.numeric(exp(predict(fit, newdata2))))
    }

    if (is.null(newdata)) {
        # Should we use the whole case-base dataset or the original one?
        # absRisk <- mean(1.0 - surv)
        # return(absRisk)
        return(NULL)
    } else {
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

        absRisk <- 1.0 - surv
        names(absRisk) <- rownames(newdata)
        return(absRisk)
    }
}





