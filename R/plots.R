#' Population Time Plot Data
#'
#' Create a data frame for population time plots to give a visual representation
#' of incidence density
#'
#' @param data a \code{data.frame} or \code{data.table} containing the source
#'   dataset.
#' @param time a character string giving the name of the time variable. See
#'   Details.
#' @param event a character string giving the name of the event variable
#'   contained in \code{data}. See Details. If \code{event} is a numeric
#'   variable, then 0 needs to represent a censored observation, 1 needs to be
#'   the event of interest. Integers 2, 3, ... and so on are treated as
#'   competing events. If event is a \code{factor} or \code{character} and
#'   \code{censored.indicator} is not specified, this function will assume the
#'   reference level is the censored indicator
#' @param censored.indicator a character string of length 1 indicating which
#'   value in \code{event} is the censored. This function will use
#'   \code{\link[stats]{relevel}} to set \code{censored.indicator} as the
#'   reference level. This argument is ignored if the \code{event} variable is a
#'   numeric
#' @param exposure a character string of length 1 giving the name of the
#'   exposure variable which must be contained in \code{data}. Default is
#'   \code{NULL}. This is used to produced exposure stratified plots. If an
#'   \code{exposure} is specified, \code{popTime} returns an object of class
#'   \code{popTimeExposure}
#' @param percentile_number Default=0.5. Give a value between 0-1. if the
#'   percentile number of available subjects at any given point is less than 10,
#'   then sample regardless of case status. Depending on distribution of
#'   survival times and events event points may not be evenly distributed with
#'   default value.
#'
#'
#' @details It is assumed that \code{data} contains the two columns
#'   corresponding to the supplied time and event variables. If either the
#'   \code{time} or \code{event} argument is missing, the function looks for
#'   columns that contain the words \code{"time"}, \code{"event"}, or
#'   \code{"status"} in them (case insensitive). The function first looks for
#'   the time variable, then it looks for the event variable. This order of
#'   operation is important if for example the time variable is named
#'   \code{"event time"} and the event variable is named \code{"event
#'   indicator"}. This function will first (automatically) find the time
#'   variable and remove this as a possibility from subsequent searches of the
#'   event variable. The following regular expressions are used for the time and
#'   event variables: \describe{ \item{time}{\code{"[\\s\\W_]+time|^time\\b"}}
#'   \item{event}{\code{"[\\s\\W_]+event|^event\\b|[\\s\\W_]+status|^status\\b"}}
#'    } This allows for \code{"time"} to be preceded or followed by one or more
#'   white space characters, one or more non-word characters or one or more
#'   underscores. For example, the following column names would be recognized by
#'   the function as the \code{"time"} variable: \code{"time of death",
#'   "death_time", "Time", "time", "diagnosis_time", "time.diag", "diag__time"}.
#'   But the following will not be recognized: \code{"diagtime","eventtime",
#'   "Timediag"}
#' @return An object of class \code{popTime} (or \code{popTimeExposure} if
#'   exposure is specified), \code{data.table} and \code{data.frame} in this
#'   order! The output of this function is to be used with the plot method for
#'   objects of class \code{popTime} or of class \code{popTimeExposure}, which
#'   will produce population time plots. This dataset augments the original data
#'   with the following columns: \describe{\item{original.time}{value of the
#'   time variable in the original dataset - the one specified by the
#'   \code{time} user argument to this function}\item{original.event}{value of
#'   the event variable in the original dataset - the one specified by the
#'   \code{event} user argument to this function}\item{time}{renames the user
#'   specified time column to time}\item{event}{renames the user specified event
#'   argument to event}}
#' @seealso \code{\link{plot.popTime}}
#' @examples
#' data("bmtcrr")
#' popTimeData <- popTime(data = bmtcrr, time = "ftime")
#' class(popTimeData)
#' popTimeData <- popTime(data = bmtcrr, time = "ftime", exposure = "D")
#' attr(popTimeData, "exposure")
#' @import data.table
#' @export
popTime <- function(data, time, event, censored.indicator,
                    exposure,percentile_number){

    #data <- DTsim;censored.indicator = NULL;exposure = "z" ; time = "time"; event = "event"
    #names(data)
    varNames <- checkArgsTimeEvent(data = data, time = time, event = event)
    #varNames <- checkArgsTimeEvent(data)
    #varNames <- checkArgsTimeEvent(data, time = time)
    ycoord <- yc <- n_available <- NULL

    DT <- data.table::as.data.table(data)
    if (missing(percentile_number)) {
        percentile_number <- 0.5
    }
    if (missing(censored.indicator)) {
        censored.indicator <- NULL
    }
    if (missing(exposure)) {

        nobs <- nrow(DT)

        # names(veteran)
        DT[, "original.time" := get(varNames$time)]
        DT[, "original.event" := get(varNames$event)]

        if (varNames$time != "time") setnames(DT,varNames$time, "time")
        if (varNames$event != "event") setnames(DT,varNames$event, "event")
# browser()
        modifiedEvent <- checkArgsEventIndicator(data = data, event = varNames$event,
                                                 censored.indicator = censored.indicator)

        DT[, event := modifiedEvent$event.numeric]
        DT[, "event status" := modifiedEvent$event.factored]
        # nLevels <- modifiedEvent$nLevels

        # people with
        # short values of t at the top
        DT[ DT[,order(time)], ycoord := (nobs:1)]

        # sample y coordinates for each event, so that we can see the incidence density
        # on population-time plots. Sampling from people who have an
        # observed time t greater than that of a fixed individual who had the event

        # need to
        # if there are only two levels, then find out how many controls
        # are left to sample from. if there are three levels,
        # check to see if there are enough 0's and 2's to sample from (this
        # implicitly assumes event=1 is the event of interest)
        # we only plot events==1 (i.e. the event of interest)
        DT[, yc := 0L]
        DT[, n_available := 0L]

        DT[event == 1, n_available := sapply(time,
                                        function(i) DT[time >= i & event != 1, .N])]

        # if the 50th percentile number of available subjects at any given
        # point is less than 10, then sample regardless of case status
        ###NEED TO MAKE THIS LESS STRINGENT##############??????
        if (DT[,quantile(n_available, probs = percentile_number)] < 15) {
            # message("Sampling from all remaining individuals under study,
            #         regardless of event status")
            DT[event == 1,
               n_available := sapply(time,
                                     function(i) DT[time >= i , .N ])]

            DT[event == 1 & n_available > 0,
               yc := sapply(time,
                            function(i)
                                sample(DT[time >= i, ycoord], 1))]

            # use original coordinate if there is no one left to sample from
            DT[event == 1 & n_available == 0, yc := ycoord]
        } else {

            # message("Sampling only from individuals who never experienced
            #         the event of interest")

            DT[event == 1 & n_available > 0,
               yc := sapply(time,
                            function(i)
                                sample(DT[time >= i & event != 1, ycoord], 1))]

            # use original coordinate if there is no one left to sample from
            DT[event == 1 & n_available == 0, yc := ycoord]

        }
        class(DT) <- c("popTime",class(DT))
        attr(DT, "exposure") <- NULL
        attr(DT, "call") <- match.call()
        return(DT)

    } else {

        DT[, "original.time" := get(varNames$time)]
        DT[, "original.event" := get(varNames$event)]

        if (varNames$time != "time") setnames(DT,varNames$time, "time")
        if (varNames$event != "event") setnames(DT,varNames$event, "event")

        l <- split(DT, DT[[exposure]])
        l <- lapply(l,
                function(i) {
                    transform(i,
                    event = checkArgsEventIndicator(data = i, event = "event",
                          censored.indicator = censored.indicator)$event.numeric,
                   `event status` = checkArgsEventIndicator(data = i, event = "event",
                              censored.indicator = censored.indicator)$event.factor
                        )
                    }
        )

        lapply(l, function(i) {
            nobs <- nrow(i)
            i[ i[,order(time)], ycoord := (nobs:1)]
        })

        # sample y coordinates for each event, so that we can see the incidence density
        # on population-time plots. Sampling from people who have an
        # observed time t greater than that of a fixed individual who had the event
        # if there are only two levels, then find out how many controls
        # are left to sample from. if there are three levels,
        # check to see if there are enough 0's and 2's to sample from (this
        # implicitly assumes event=1 is the event of interest)
        # we only plot events==1 (i.e. the event of interest)

        lapply(l, function(K) {
            K[, yc := 0L]
            K[, n_available := 0L]

            K[event == 1, n_available := sapply(time,
                                                 function(i) K[time >= i & event != 1, .N])]

            # if the 50th percentile number of available subjects at any given
            # point is less than 10, then sample regardless of case status
            if (K[,quantile(n_available, probs = percentile_number)] < 10) {
                # message("Sampling from all remaining individuals under study,
                #     regardless of event status")
                K[event == 1,
                   n_available := sapply(time,
                                         function(i) K[time >= i , .N ])]

                K[event == 1 & n_available > 0,
                   yc := sapply(time,
                                function(i)
                                    sample(K[time >= i, ycoord], 1))]

                # use original coordinate if there is no one left to sample from
                K[event == 1 & n_available == 0, yc := ycoord]
            } else {

                # message("Sampling only from individuals who never experienced
                #     the event of interest")

                K[event == 1 & n_available > 0,
                   yc := sapply(time,
                                function(i)
                                    sample(K[time >= i & event != 1, ycoord], 1))]

                # use original coordinate if there is no one left to sample from
                K[event == 1 & n_available == 0, yc := ycoord]

            }
        }
        )

        lk <- data.table::rbindlist(l)
        attr(lk, "exposure") <- exposure
        class(lk) <- c("popTime", class(lk))
        attr(lk, "call") <- match.call()
        return(lk)

    }

}

# taken verbatim from cowplot::theme_cowplot()
theme_cb <- function(font_size = 14, font_family = "", line_size = 0.5,
          rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14) {
    half_line <- font_size/2
    small_size <- rel_small * font_size
    theme_grey(base_size = font_size, base_family = font_family) %+replace%
        theme(line = element_line(color = "black", size = line_size,
                                  linetype = 1, lineend = "butt"),
              rect = element_rect(fill = NA,
                                  color = NA, size = line_size, linetype = 1),
              text = element_text(family = font_family,
                                  face = "plain", color = "black", size = font_size,
                                  hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                  margin = margin(), debug = FALSE),
              axis.line = element_line(color = "black",
                                       size = line_size, lineend = "square"),
              axis.line.x = NULL,
              axis.line.y = NULL,
              axis.text = element_text(color = "black",
                                       size = small_size),
              axis.text.x = element_text(margin = margin(t = small_size/4),
                                         vjust = 1),
              axis.text.x.top = element_text(margin = margin(b = small_size/4),
                                             vjust = 0),
              axis.text.y = element_text(margin = margin(r = small_size/4),
                                         hjust = 1),
              axis.text.y.right = element_text(margin = margin(l = small_size/4),
                                               hjust = 0),
              axis.ticks = element_line(color = "black",
                                        size = line_size),
              axis.ticks.length = unit(half_line/2,
                                       "pt"),
              axis.title.x = element_text(margin = margin(t = half_line/2),
                                          vjust = 1),
              axis.title.x.top = element_text(margin = margin(b = half_line/2),
                                              vjust = 0),
              axis.title.y = element_text(angle = 90,
                                          margin = margin(r = half_line/2), vjust = 1),
              axis.title.y.right = element_text(angle = -90, margin = margin(l = half_line/2),
                                                vjust = 0),
              legend.background = element_blank(),
              legend.spacing = unit(font_size, "pt"),
              legend.spacing.x = NULL,
              legend.spacing.y = NULL,
              legend.margin = margin(0,
                                     0, 0, 0),
              legend.key = element_blank(),
              legend.key.size = unit(1.1 *
                                         font_size, "pt"),
              legend.key.height = NULL,
              legend.key.width = NULL,
              legend.text = element_text(size = rel(rel_small)),
              legend.text.align = NULL,
              legend.title = element_text(hjust = 0),
              legend.title.align = NULL,
              legend.position = "right",
              legend.direction = NULL,
              legend.justification = c("left",
                                       "center"),
              legend.box = NULL, legend.box.margin = margin(0,
                                                            0, 0, 0),
              legend.box.background = element_blank(),
              legend.box.spacing = unit(font_size, "pt"),
              panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              panel.grid.major = NULL,
              panel.grid.minor = NULL,
              panel.grid.major.x = NULL,
              panel.grid.major.y = NULL,
              panel.grid.minor.x = NULL,
              panel.grid.minor.y = NULL,
              panel.spacing = unit(half_line, "pt"),
              panel.spacing.x = NULL,
              panel.spacing.y = NULL,
              panel.ontop = FALSE,
              strip.background = element_rect(fill = "grey80"),
              strip.text = element_text(size = rel(rel_small),
                                        margin = margin(half_line/2, half_line/2, half_line/2,
                                                        half_line/2)),
              strip.text.x = NULL,
              strip.text.y = element_text(angle = -90),
              strip.placement = "inside",
              strip.placement.x = NULL,
              strip.placement.y = NULL,
              strip.switch.pad.grid = unit(half_line/2,
                                           "pt"),
              strip.switch.pad.wrap = unit(half_line/2,
                                           "pt"),
              plot.background = element_blank(),
              plot.title = element_text(face = "bold",
                                        size = rel(rel_large), hjust = 0, vjust = 1,
                                        margin = margin(b = half_line)),
              plot.subtitle = element_text(size = rel(rel_small),
                                           hjust = 0, vjust = 1, margin = margin(b = half_line)),
              plot.caption = element_text(size = rel(rel_tiny),
                                          hjust = 1, vjust = 1, margin = margin(t = half_line)),
              plot.tag = element_text(face = "bold", hjust = 0,
                                      vjust = 0.7), plot.tag.position = c(0, 1),
              plot.margin = margin(half_line,
                                   half_line, half_line, half_line), complete = TRUE)
}



# taken verbatim from cowplot::panel_border
panelBorder <- function(color = "grey85", size = 1, linetype = 1, remove = FALSE,
                        colour) {
    if (!missing(colour)) {
        color <- colour
    }
    if (remove) {
        return(theme(panel.border = element_blank()))
    }
    theme(panel.border = element_rect(color = color, fill = NA,
                                      linetype = linetype, size = size))
}



#' @rdname plot.singleEventCB
plotHazardRatio <- function(x, newdata, newdata2, ci, ci.lvl, ci.col,
                            rug, xvar, ...) {

    tt <- terms(x)
    Terms <- delete.response(tt)
    beta2 <- coef(x)

    gradient <- hrJacobian(object = x, newdata = newdata,
                           newdata2 = newdata2, term = Terms)

    log_hazard_ratio <- gradient %*% beta2

    if (is.null(xvar)){
        xvar <- x[["timeVar"]]
    } else {
        if (length(xvar) > 1) warning("more than one xvar supplied. Only plotting hazard ratio for first element.")
        xvar <- xvar[1]
    }

    if (xvar %ni% colnames(newdata))
        stop(sprintf("%s column (which you supplied to 'xvar' argument) not found in newdata",xvar)) else
            xvar_values <- newdata[[xvar]]

    # browser()
    # sorting indices for ploting
    i.backw <- order(xvar_values, decreasing = TRUE)
    i.forw <- order(xvar_values)

    # plot CI as polygon shade - if 'se = TRUE' (default)
    if (ci) {

        v2 <- vcov(x)
        SE_log_hazard_ratio <- sqrt(diag(gradient %*% tcrossprod(v2, gradient)))

        hazard_ratio_lower <- exp(qnorm(p = (1 - ci.lvl) / 2, mean = log_hazard_ratio, sd = SE_log_hazard_ratio))
        hazard_ratio_upper <- exp(qnorm(p = 1 - (1 - ci.lvl) / 2, mean = log_hazard_ratio, sd = SE_log_hazard_ratio))
        x.poly <- c(xvar_values[i.forw] , xvar_values[i.backw])
        y.poly <- c(hazard_ratio_lower[i.forw] , hazard_ratio_upper[i.backw])

        other_plot_args <- list(...)

        line_args <- grep("^lwd$|^lty$|^col$", names(other_plot_args))

        if (length(line_args) == 0)
            line_args <- list(NULL) else
                line_args <- other_plot_args[line_args]

        do.call("plot", utils::modifyList(
            list(x = range(x.poly),
                 y = range(y.poly),
                 type = "n",
                 ylab = "hazard ratio",
                 xlab = xvar),
            other_plot_args
        ))

        polygon(x.poly , y.poly , col = ci.col , border = NA)

        do.call("lines", utils::modifyList(
            list(x = xvar_values,
                 y = exp(log_hazard_ratio),
                 lwd = 2,
                 lty = 1,
                 col = "black"),
            line_args
        ))

        # browser()

        results <- transform(newdata,
                             log_hazard_ratio = log_hazard_ratio,
                             standarderror = SE_log_hazard_ratio,
                             hazard_ratio = exp(log_hazard_ratio),
                             lowerbound = hazard_ratio_lower,
                             upperbound = hazard_ratio_upper)

    } else {
        do.call("plot", utils::modifyList(
            list(x = xvar_values, y = exp(log_hazard_ratio), lwd = 2, lty = 1, type = "l",
                 ylab = "hazard ratio", xlab = xvar),
            list(...)
        ))

        results <- transform(newdata,
                             log_hazard_ratio = log_hazard_ratio,
                             hazard_ratio = exp(log_hazard_ratio))
    }

    if (rug) {
        events <- x[["originalData"]][[x[["eventVar"]]]]
        rug(x[["originalData"]][which(events==1),,drop=F][[ x[["timeVar"]]  ]],
            quiet = TRUE) # Silence warnings about clipped values
    }

    invisible(results)

}
