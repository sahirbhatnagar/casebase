#' @rdname plot.singleEventCB
#' @param newdata2 `data.frame` for exposed group. calculated and passed
#'   internally to `plotHazardRatio` function
#' @importFrom graphics polygon lines arrows points
#' @importFrom stats terms delete.response vcov qnorm
#' @importFrom utils modifyList
plotHazardRatio <- function(x, newdata, newdata2, ci, ci.lvl, ci.col,
                            rug, xvar, ...) {
    tt <- stats::terms(x)
    Terms <- stats::delete.response(tt)
    beta2 <- coef(x)
    gradient <- hrJacobian(
        object = x, newdata = newdata,
        newdata2 = newdata2, term = Terms
    )

    log_hazard_ratio <- gradient %*% beta2

    if (is.null(xvar)) {
        xvar <- x[["timeVar"]]
    } else {
        if (length(xvar) > 1) {
            warning(paste("more than one xvar supplied. Only plotting hazard",
                          "ratio for first element."))
            }
        xvar <- xvar[1]
    }

    if (xvar %ni% colnames(newdata)) {
        stop(sprintf(paste("%s column (which you supplied to 'xvar' argument)",
                           "not found in newdata"), xvar))
    } else {
        xvar_values <- newdata[[xvar]]
    }

    # sorting indices for ploting
    i.backw <- order(xvar_values, decreasing = TRUE)
    i.forw <- order(xvar_values)

    other_plot_args <- list(...)
    line_args <- grep("^lwd$|^lty$|^col$|^pch$|^cex$", names(other_plot_args))

    if (length(line_args) == 0) {
        line_args <- list(NULL)
    } else {
        line_args <- other_plot_args[line_args]
    }

    # plot CI as polygon shade - if 'se = TRUE' (default)
    if (ci) {
        v2 <- stats::vcov(x)
        SE_log_hazard_ratio <- sqrt(diag(gradient %*% tcrossprod(v2, gradient)))

        hazard_ratio_lower <- exp(stats::qnorm(p = (1 - ci.lvl) / 2,
                                               mean = log_hazard_ratio,
                                               sd = SE_log_hazard_ratio))
        hazard_ratio_upper <- exp(stats::qnorm(p = 1 - (1 - ci.lvl) / 2,
                                               mean = log_hazard_ratio,
                                               sd = SE_log_hazard_ratio))
        x.poly <- c(xvar_values[i.forw], xvar_values[i.backw])
        y.poly <- c(hazard_ratio_lower[i.forw], hazard_ratio_upper[i.backw])

        do.call("plot", utils::modifyList(
            list(
                x = range(x.poly),
                y = range(y.poly) * c(0.99, 1.01),
                type = "n",
                ylab = "hazard ratio",
                xlab = xvar
            ),
            other_plot_args
        ))

        if (length(unique(x.poly)) == 1) {
            graphics::arrows(x0 = x.poly[1], y0 = y.poly[1], y1 = y.poly[2],
                             angle = 90, length = 0.5, code = 3)
            do.call("points", utils::modifyList(
                list(
                    x = xvar_values[i.forw],
                    y = exp(log_hazard_ratio[i.forw]),
                    pch = 19,
                    col = "black",
                    cex = 2
                ),
                line_args
            ))
        } else {
            graphics::polygon(x.poly, y.poly, col = ci.col, border = NA)
            do.call("lines", utils::modifyList(
                list(
                    x = xvar_values[i.forw],
                    y = exp(log_hazard_ratio[i.forw]),
                    lwd = 2,
                    lty = 1,
                    col = "black"
                ),
                line_args
            ))
        }

        results <- transform(newdata,
                             log_hazard_ratio = log_hazard_ratio,
                             standarderror = SE_log_hazard_ratio,
                             hazard_ratio = exp(log_hazard_ratio),
                             lowerbound = hazard_ratio_lower,
                             upperbound = hazard_ratio_upper
        )

    } else {

        if (length(xvar_values) == 1) {
            do.call("plot", utils::modifyList(
                list(
                    x = xvar_values, y = exp(log_hazard_ratio), lwd = 2,
                    lty = 1, pch = 19, cex = 2, ylab = "hazard ratio",
                    xlab = xvar
                ),
                other_plot_args
            ))

        } else {

            do.call("plot", utils::modifyList(
                list(
                    x = xvar_values, y = exp(log_hazard_ratio), lwd = 2,
                    lty = 1, type = "l", ylab = "hazard ratio", xlab = xvar
                ),
                other_plot_args
            ))
        }

        results <- transform(newdata,
                             log_hazard_ratio = log_hazard_ratio,
                             hazard_ratio = exp(log_hazard_ratio)
        )
    }

    if (rug) {
        events <- x[["originalData"]][[x[["eventVar"]]]]
        rug(x[["originalData"]][which(events == 1), , drop = FALSE][[xvar]],
            quiet = TRUE
        ) # Silence warnings about clipped values
    }

    invisible(results)
}
