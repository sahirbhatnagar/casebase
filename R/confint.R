#' Compute confidence intervals for risks
#'
#' This function uses parametric bootstrap to compute confidence intervals for
#' the risk estimates. Since it relies on MLE theory for the validity of these
#' intervals, this function only works for \code{fit_obj} that was fitted using
#' \code{family = "glm"} (i.e. the default).
#'
#' If the package \code{progress} is available, the function also reports on
#' progress of the sampling (which can take some time if there are many
#' covariate profiles and/or time points).
#'
#' @param object Output of function \code{absoluteRisk}.
#' @param parm Output of function \code{fitSmoothHazard} that was used to
#'   compute \code{object}.
#' @param level The confidence level required.
#' @param nboot The number of bootstrap samples to use.
#' @param ... Additional arguments for methods.
#' @return If there is only one covariate profile, the function returns a matrix
#'   with the time points, the risk estimates, and the confidence intervals. If
#'   there are more than one covariate profile, the function returns a list with
#'   three components.
#' @importFrom stats confint
#' @export
confint.absRiskCB <- function(object, parm, level = 0.95,
                              nboot = 500, ...) {
    # Relies on MLE theory for glm
    if (!inherits(parm, "glm")) {
        stop(paste("fit_obj is of class", class(parm)[1],
                   "\naconfint should be used with a fit_obj of class glm"),
             call. = TRUE)
    }
    # If available, add progress bar
    prog_bar <- requireNamespace("progress", quietly = TRUE)
    format_pb <- "Sampling [:bar] :current/:total (:percent)"
    if (prog_bar) {
        pb <- progress::progress_bar$new(format = format_pb,
                                         total = nboot)
    } else message("Waiting for profiling to be done...")

    # Extract coef + vcov matrix
    mu <- coef(parm)
    sigma <- vcov(parm)

    # Sample from multivariate normal
    betas_boot <- generate_mvnorm(nboot, mu, sigma)
    # Extract required data
    newdata <- attr(object, "newdata")
    times <- object[,"time"]
    # Create placeholder for replacing coeffs
    fit_obj_tmp <- parm

    # For each bootstrap sample, rerun absRisk computation
    abs_boot <- lapply(seq_len(nboot), function(b) {
        # If progress bar is available, tick
        if (prog_bar) pb$tick()
        fit_obj_tmp$coefficients <- betas_boot[b, ]
        abs_tmp <- absoluteRisk(fit_obj_tmp, time = times,
                                newdata = newdata)
        return(abs_tmp[,-1])
    })

    # Clean output
    if (nrow(newdata) > 1) {
        # Compute quantiles
        quant_tmp <- lapply(seq_len(nrow(newdata)), function(i) {
            mat <- sapply(abs_boot, function(elem) elem[,i])
            extract_quantiles(mat, level)
            })
        lower_ci <- sapply(quant_tmp, function(elem) elem[1,])
        lower_ci <- cbind(times, lower_ci)
        colnames(lower_ci)[1] <- "time"
        upper_ci <- sapply(quant_tmp, function(elem) elem[2,])
        upper_ci <- cbind(times, upper_ci)
        colnames(upper_ci)[1] <- "time"
        lower_name <- create_perc(0.5 - 0.5*level)
        upper_name <- create_perc(0.5 + 0.5*level)

        out <- list(object, lower_ci, upper_ci)
        names(out) <- c("Estimate", lower_name, upper_name)

    } else {
        # If only one row, add quantiles next to estimates
        abs_boot <- simplify2array(abs_boot)
        out <- extract_quantiles(abs_boot, level)
        out <- cbind(object, t(out))
    }

    return(out)
}

# Helper Functions----
#' @importFrom stats rnorm
generate_mvnorm <- function(n, mu, sigma) {
    p <- ncol(sigma)
    # Compute sqrt of sigma
    chol_decomp <- chol(sigma, pivot = TRUE)
    R <- chol_decomp[, order(attr(chol_decomp, "pivot"))]
    # Generate data with right covariance
    data <- matrix(rnorm(n*p), ncol = p) %*% R
    # Add back mean
    data <- sweep(data, 2, mu, "+")
    return(data)
}

extract_quantiles <- function(mat, level) {
    apply(mat, 1, function(row) {
        quantile(row, probs = c(0.5 - 0.5*level,
                                0.5 + 0.5*level))
    })
}

create_perc <- function(x) {
    paste0(format(100*x, scientific = FALSE, trim = TRUE,
                  digits = max(2L, getOption("digits"))),
           "%")
}
