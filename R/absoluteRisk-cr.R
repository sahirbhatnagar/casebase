#' @rdname absoluteRisk
#' @importFrom data.table := data.table
#' @export
absoluteRisk.CompRisk <- function(object, time, newdata,
                                  method = c("numerical", "montecarlo"),
                                  nsamp = 100, onlyMain = TRUE,
                                  type = c("CI", "survival"),
                                  addZero = TRUE) {
    method <- match.arg(method)
    type <- match.arg(type)

    if (missing(newdata)) {
        # Should we use the whole case-base dataset or the original one?
        if (is.null(object@originalData)) {
            stop(paste("Cannot estimate the mean absolute risk without",
                       "the original data. See documentation."),
                 call. = FALSE)
        }
        newdata <- object@originalData
        unselectTime <- (names(newdata) != object@timeVar)
        newdata <- subset(newdata, select = unselectTime)
    }
    typeEvents <- object@typeEvents
    ###################################################
    # In competing risks, we can get a cumulative
    # incidence function using a nested double integral
    # subdensity f_j = lambda_j * Survival
    # F_j = P(T <= t, J = j : covariates) = int_0^t f_j
    ###################################################
    time_ordered <- unique(c(0, sort(time)))
    # Create array to store output
    if (onlyMain) {
        output <- matrix(NA, ncol = nrow(newdata) + 1,
                         nrow = length(time_ordered))
        output[, 1] <- time_ordered
        colnames(output) <- c("time", rep("", nrow(newdata)))
        rownames(output) <- rep("", length(time_ordered))
        output[1, -1] <- 0
    } else {
        J <- length(typeEvents) - 1
        output <- array(NA, dim = c(length(time_ordered), nrow(newdata) + 1, J))
        output[, 1, ] <- time_ordered
        dimnames(output) <- list(rep("", length(time_ordered)),
                                 c("time", rep("", nrow(newdata))),
                                 paste("event", typeEvents[-1], sep = "="))
        output[1, -1, ] <- 0
    }

    # Compute subdensities
    if (method == "numerical") {
        # Compute points at which we evaluate integral
        knots <- seq(0, max(time_ordered),
                     length.out = (length(time_ordered) - 1) * nsamp)
        knots <- unique(sort(c(knots, time_ordered)))

        for (j in seq_len(nrow(newdata))) {
            # Extract current obs
            current_obs <- newdata[j, , drop = FALSE]
            # Create data.table for prediction
            newdata2 <- data.table::data.table(current_obs)
            newdata2 <- newdata2[rep(1, length(knots))]
            newdata2[, object@timeVar := knots]
            newdata2[, "offset" := 0]
            # Compute all values for all hazards
            lambdas <- exp(predict_CompRisk(object, newdata2))
            lambdas[which(lambdas %in% c(Inf, -Inf))] <- 0
            OverallLambda <- rowSums(lambdas)
            survFunction <- exp(-trap_int(knots, OverallLambda))
            if (onlyMain) {
                # Only compute first subdensity
                subdensity <- lambdas[, 1] * survFunction
                pred <- trap_int(knots, subdensity)[knots %in% c(0, time)]
                output[, j + 1] <- pred
            } else {
                subdensity <- lambdas * drop(survFunction)
                pred <- trap_int(knots, subdensity)[knots %in% c(0, time)]
                output[, j + 1, ] <- pred
            }
        }
    } else {
        # Sample points at which we evaluate function
        knots <- runif(n = (length(time_ordered) - 1) * nsamp^2,
                       min = 0, max = max(time_ordered))
        knots2 <- runif(n = (length(time_ordered) - 1) * nsamp,
                        min = 0, max = max(time_ordered))
        knots2 <- sort(knots2)
        for (j in seq_len(nrow(newdata))) {
            current_obs <- newdata[j, , drop = FALSE]
            # Create data.table for prediction
            newdata2 <- data.table(current_obs)
            newdata2 <- newdata2[rep(1, length(knots))]
            newdata2[, object@timeVar := knots]
            newdata2[, "offset" := 0]
            # Compute all values for all hazards
            lambdas <- exp(predict_CompRisk(object, newdata2))
            lambdas[which(lambdas %in% c(Inf, -Inf))] <- 0
            OverallLambda <- rowSums(lambdas)
            mean_values <- sapply(split(OverallLambda,
                                        cut(knots, breaks = knots2)),
                                  mean, na.rm = TRUE)
            mean_values[is.na(mean_values)] <- 0
            survFunction <- exp(-cumsum(c(0, mean_values * diff(knots2))))
            # Second integral---Create data.table for prediction
            newdata2 <- data.table(current_obs)
            newdata2 <- newdata2[rep(1, length(knots2))]
            newdata2[, object@timeVar := knots2]
            newdata2[, "offset" := 0]
            # Compute all values for all hazards
            lambdas2 <- exp(predict_CompRisk(object, newdata2))
            lambdas2[which(lambdas2 %in% c(Inf, -Inf))] <- 0
            if (onlyMain) {
                # Only compute first subdensity
                subdensity <- lambdas2[, 1] * survFunction
                mean_values2 <- sapply(split(subdensity,
                                             cut(knots2,
                                                 breaks = time_ordered)),
                                       mean, na.rm = TRUE)
                pred <- cumsum(c(0, mean_values2 * diff(time_ordered)))
                output[, j + 1] <- pred
            } else {
                subdensity <- lambdas2 * drop(survFunction)
                mean_values2 <- sapply(split(seq_len(nrow(subdensity)),
                                             cut(knots2,
                                                 breaks = time_ordered)),
                                       function(ind) {
                                           colMeans(subdensity[ind, ,
                                                               drop = FALSE],
                                                    na.rm = TRUE)
                                           })
                pred <- apply(mean_values2, 1, function(row) {
                    cumsum(c(0, row * diff(time_ordered)))
                    })
                output[, j + 1, ] <- pred
            }
        }
    }

    if (onlyMain) {
        # Switch to survival scale?
        if (type == "survival") {
            output[, -1] <- 1 - output[, -1, drop = FALSE]
        }
        # Add time = 0?
        if (!addZero) output <- output[-1, , drop = FALSE]
    } else {
        # Switch to survival scale?
        if (type == "survival") {
            output[, -1, ] <- 1 - output[, -1, ]
        }
        # Add time = 0?
        if (!addZero) output <- output[-1, , , drop = FALSE]

    }
    # Sometimes montecarlo integration gives nonsensical probability estimates
    if (method == "montecarlo" && (any(output < 0) | any(output > 1))) {
        warning(paste("Some probabilities are out of range. Consider",
                      "increasing nsamp or using numerical integration"),
                call. = FALSE)
    }
    attr(output, "type") <- type
    return(output)

}

# #' @rdname absoluteRisk
# #' @export
absoluteRisk.CompRiskGlmnet <- function(object, time, newdata,
                                        method = c("numerical", "montecarlo"),
                                        nsamp = 100, onlyMain = TRUE,
                                        s = c("lambda.1se", "lambda.min"),
                                        type = c("CI", "survival"),
                                        addZero = TRUE, ...) {
    # The current implementation doesn't work
    stop(paste("object is of class", class(object),
               "\nabsoluteRisk should be used with an object of class glm,",
               "cv.glmnet, gbm, or CompRisk"),
         call. = TRUE)
}

#' @importFrom stats coef
predict_CompRisk <- function(object, newdata = NULL) {
    ttob <- terms(object)
    X <- model.matrix(delete.response(ttob), newdata,
                      contrasts = if (length(object@contrasts)) {
                          object@contrasts
                          } else NULL,
                      xlev = object@xlevels)
    coeffs <- matrix(coef(object), nrow = ncol(X),
                     byrow = TRUE)
    preds <- X %*% coeffs
    colnames(preds) <- paste0("log(mu[,",
                              seq(2, length(object@typeEvents)),
                              "]/mu[,1])")
    return(preds)
}
