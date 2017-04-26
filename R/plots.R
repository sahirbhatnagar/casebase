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
#'   event variable.
#'   The following regular expressions are used for the time and event
#'   variables: \describe{ \item{time}{\code{"[\\s\\W_]+time|^time\\b"}}
#'   \item{event}{\code{"[\\s\\W_]+event|^event\\b|[\\s\\W_]+status|^status\\b"}}
#'    }
#'   This allows for \code{"time"} to be preceded or followed by one or more
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
#'   will produce population time plots
#' @seealso \code{\link{plot.popTime}}, \code{\link{plot.popTimeExposure}}
#'
#' @import data.table
#' @export
popTime <- function(data, time, event, censored.indicator,
                    exposure){

    #data <- DTsim;censored.indicator = NULL;exposure = "z" ; time = "time"; event = "event"
    #names(data)
    varNames <- checkArgsTimeEvent(data = data, time = time, event = event)
    #varNames <- checkArgsTimeEvent(data)
    #varNames <- checkArgsTimeEvent(data, time = time)
    ycoord <- yc <- n_available <- NULL

    DT <- data.table::as.data.table(data)
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
        if (DT[,quantile(n_available, probs = 0.5)] < 10) {
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
                    event = checkArgsEventIndicator(data = i, event = varNames$event,
                          censored.indicator = censored.indicator)$event.numeric,
                   `event status` = checkArgsEventIndicator(data = i, event = varNames$event,
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
            if (K[,quantile(n_available, probs = 0.5)] < 10) {
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
        lkj <- list(data = lk, exposure = exposure)
        class(lkj) <- c("popTimeExposure", class(lkj))
        attr(lkj, "call") <- match.call()
        return(lkj)

    }

}
