#' @rdname absoluteRisk
#' @export
absoluteRisk.CompRisk <- function(object, time, newdata, method = c("numerical", "montecarlo"),
                                  nsamp = 100, onlyMain = TRUE, ...) {
    method <- match.arg(method)

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
    }
    timeVar <- object@timeVar
    typeEvents <- object@typeEvents
    ###################################################
    # In competing risks, we can get a cumulative
    # incidence function using a nested double integral
    # f_j = lambda_j * Survival
    # F_j = P(T <= t, J = j : covariates) = int_0^t f_j
    ###################################################

    lambda_vec <- function(x, object, newdata) {
        # Note: the offset should be set to zero when estimating the hazard.
        newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                               row.names = as.character(1:length(x)))
        newdata2[timeVar] <- x
        # predictvglm doesn't like offset = 0
        withCallingHandlers(pred <- VGAM::predictvglm(object, newdata2),
                            warning = handler_offset)
        return(exp(pred))
    }

    # Compute subdensities
    if (method == "numerical") {
        overallLambda <- function(x, object, newdata, lambda_vec) {
            return(as.numeric(rowSums(lambda_vec(x, object, newdata))))
        }
        overallSurv <- function(x, object, newdata, lambda_vec,
                                method, nsamp) {
            fn <- function(t) {
                exp(-integrate(overallLambda, lower = 0, upper = t,
                               object = object, newdata = newdata,
                               lambda_vec = lambda_vec,
                               subdivisions = nsamp)$value)
            }
            # Vectorize
            return(sapply(x, fn))
        }
        J <- length(typeEvents) - 1
        # 2. Compute individual subdensities f_j
        subdensities <- vector("list", length = J)
        subdensity_template <- function(x, object, newdata, index) {
            newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                                   row.names = as.character(1:length(x)))
            newdata2[timeVar] <- x
            lambdas <- lambda_vec(x, object, newdata2)
            lambdas[,index] * overallSurv(x, object = object, newdata2, lambda_vec, method, nsamp)
        }
        for (j in 1:J) {
            subdensities[[j]] <- partialize(subdensity_template, index = j)
        }
    } else subdensities <- NULL

    return(estimate_risk_cr(object, time, newdata, method,
                            nsamp, onlyMain, subdensities,
                            lambda_vec, timeVar, typeEvents))

}

#' @rdname absoluteRisk
#' @export
absoluteRisk.CompRiskGlmnet <- function(object, time, newdata, method = c("numerical", "montecarlo"),
                                        nsamp = 100, onlyMain = TRUE, s = c("lambda.1se","lambda.min"), ...) {
    warning("There is currently no guarantee that the output is correct.")
    method <- match.arg(method)
    if (is.numeric(s))
        s <- s[1]
    else if (is.character(s)) {
        s <- match.arg(s)
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
    }
    timeVar <- object$timeVar
    typeEvents <- object$typeEvents
    ###################################################
    # In competing risks, we can get a cumulative
    # incidence function using a nested double integral
    # f_j = lambda_j * Survival
    # F_j = P(T <= t, J = j : covariates) = int_0^t f_j
    ###################################################
    if (is.null(object$matrix.fit)) {
        lambda_vec <- function(x, object, newdata) {
            newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                                   row.names = as.character(1:length(x)))
            newdata2[timeVar] <- x
            formula_pred <- update(formula(delete.response(terms(object$formula))), ~ . -1)
            newdata_matrix <- model.matrix(formula_pred, newdata2)
            pred <- predict(object, newdata_matrix, s, type = "response")
            pred <- exp(pred[,seq(2, ncol(pred))] - pred[,1])
            return(pred)
        }
    } else {
        lambda_vec <- function(x, object, newdata) {
            newdata_matrix <- newdata[,colnames(newdata) != timeVar, drop = FALSE]
            newdata_matrix <- as.matrix(cbind(as.data.frame(newdata_matrix),
                                              model.matrix(update(object$formula_time, ~ . -1),
                                                           setNames(data.frame(x), timeVar)),
                                              row.names = NULL))
            pred <- predict(object, newdata_matrix, s, type = "response")
            pred <- exp(pred[,seq(2, ncol(pred))] - pred[,1])
            return(pred)
        }
    }

    # Compute subdensities
    if (method == "numerical") {
        overallLambda <- function(x, object, newdata, lambda_vec) {
            return(as.numeric(rowSums(lambda_vec(x, object, newdata))))
        }
        overallSurv <- function(x, object, newdata, lambda_vec,
                                method, nsamp) {
            fn <- function(t) {
                exp(-integrate(overallLambda, lower = 0, upper = t,
                               object = object, newdata = newdata,
                               lambda_vec = lambda_vec,
                               subdivisions = nsamp)$value)
            }
            # Vectorize
            return(sapply(x, fn))
        }
        J <- length(typeEvents) - 1
        # 2. Compute individual subdensities f_j
        subdensities <- vector("list", length = J)
        if (is.null(object$matrix.fit)) {
            subdensity_template <- function(x, object, newdata, index) {
                newdata2 <- data.frame(newdata, offset = rep_len(0, length(x)),
                                        row.names = as.character(1:length(x)))
                newdata2[timeVar] <- x
                lambdas <- lambda_vec(x, object, newdata2)
                lambdas[,index] * overallSurv(x, object = object, newdata2, lambda_vec, method, nsamp)
            }
        } else {
            subdensity_template <- function(x, object, newdata, index) {
                newdata_matrix <- newdata[,colnames(newdata) != timeVar, drop = FALSE]
                newdata_matrix <- as.matrix(cbind(as.data.frame(newdata_matrix),
                                                  model.matrix(update(object$formula_time, ~ . -1),
                                                               setNames(data.frame(x), timeVar)),
                                                  row.names = NULL))
                lambdas <- lambda_vec(x, object, newdata_matrix)
                lambdas[,index] * overallSurv(x, object = object, newdata_matrix, lambda_vec, method, nsamp)
            }
        }

        for (j in 1:J) {
            subdensities[[j]] <- partialize(subdensity_template, index = j)
        }
    } else subdensities <- NULL

    return(estimate_risk_cr(object, time, newdata, method,
                            nsamp, onlyMain, subdensities,
                            lambda_vec, timeVar, typeEvents))

}

predict.CompRiskGlmnet <- function(object, ...) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Pkg glmnet needed for this function to work. Please install it.",
             call. = FALSE)
    } else drop(glmnet::predict.cv.glmnet(object, ...))
}

estimate_risk_cr <- function(object, time, newdata, method,
                             nsamp, onlyMain, subdensities,
                             lambda_vec, timeVar, typeEvents) {
    # Create output matrix
    J <- length(typeEvents) - 1
    time_ordered <- unique(c(0, sort(time)))
    output <- array(NA, dim = c(length(time_ordered), nrow(newdata) + 1, J))
    output[,1,] <- time_ordered
    dimnames(output) <- list(rep("", length(time_ordered)),
                             c("time", rep("", nrow(newdata))),
                             paste("event", typeEvents[-1], sep = "="))
    output[1,-1,] <- 0

    # 1. Compute overall survival
    overallLambda <- function(x, object, newdata, lambda_vec) {
        return(as.numeric(rowSums(lambda_vec(x, object, newdata))))
    }
    overallSurv <- function(x, object, newdata, lambda_vec,
                            method, nsamp) {
        # Need to integrate later, so need to ability to take a vector x
        if (method == "numerical") {
            fn <- function(t) {
                exp(-integrate(overallLambda, lower = 0, upper = t,
                               object = object, newdata = newdata,
                               lambda_vec = lambda_vec,
                               subdivisions = nsamp)$value)
            }
        } else if (method == "montecarlo"){
            fn <- function(t) {
                exp(-integrate_mc(overallLambda, lower = 0, upper = t,
                                  object = object, newdata = newdata,
                                  lambda_vec = lambda_vec,
                                  subdivisions = nsamp))
            }
        } else {
            stop("Unrecognised integration method")
        }
        # Vectorize
        return(sapply(x, fn))
    }

    if (length(time_ordered) > 1) {
        if (method == "numerical") {
            # 3. Compute cumulative incidence functions F_j
            for (i in 1:nrow(newdata)) {
                for (j in 1:J) {
                    for (k in 2:length(time_ordered)) {
                        output[k,i + 1,j] <- integrate(subdensities[[j]], lower = time_ordered[k - 1],
                                                       upper = time_ordered[k],
                                                       object = object, newdata = newdata[i,,drop = FALSE],
                                                       subdivisions = nsamp)$value
                    }
                    output[,i + 1,j] <- cumsum(output[,i + 1,j])
                }
            }
        }

        if (method == "montecarlo") {
            sampledPoints <- runif(nsamp)
            for (i in 1:nrow(newdata)) {
                for (k in 2:length(time_ordered)) {
                    # Integrate all subdensities at the same time
                    x_vect <- sampledPoints * (time_ordered[k] - time_ordered[k - 1]) + time_ordered[k - 1]
                    surv <- overallSurv(x_vect, object, newdata[i,,drop = FALSE], lambda_vec, method, nsamp)
                    lambdas <- lambda_vec(x_vect, object, newdata[i,,drop = FALSE])
                    output[k, i + 1, ] <- (time_ordered[k] - time_ordered[k - 1]) * colMeans(surv * lambdas)
                }
                # if k==1, there was only one time point and we don't need to sum the contributions
                if (k != 1) output[,i + 1,] <- apply(output[,i + 1,], 2, cumsum)
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
    # Sometimes montecarlo integration gives nonsensical probability estimates
    if (method == "montecarlo" && (any(output < 0) | any(output > 1))) {
        warning("Some probabilities are out of range. Consider increasing nsamp or using numerical integration", call. = FALSE)
    }
    # If there is only one time point, we should drop a dimension and return a matrix
    if (onlyMain) return(output[,,1, drop = TRUE]) else return(output)
}
