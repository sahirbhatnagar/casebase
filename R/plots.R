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
#'
#'   The following regular expressions are used for the time and event
#'   variables: \describe{ \item{time}{\code{"[\\s\\W_]+time|^time\\b"}}
#'   \item{event}{\code{"[\\s\\W_]+event|^event\\b|[\\s\\W_]+status|^status\\b"}}
#'    }
#'
#'   This allows for \code{"time"} to be preceded or followed by one or more
#'   white space characters, one or more non-word characters or one or more
#'   underscores. For example, the following column names would be recognized by
#'   the function as the \code{"time"} variable: \code{"time of death",
#'   "death_time", "Time", "time", "diagnosis_time", "time.diag", "diag__time"}.
#'   But the following will not be recognized: \code{"diagtime","eventtime",
#'   "Timediag"}
#'
#' @return an object of class \code{popTime} (or \code{popTimeExposure} if
#'   exposure is specified), \code{data.table} and \code{data.frame} in this
#'   order! The output of this function is to be used with the plot method for
#'   objects of class \code{popTime} or of class \code{popTimeExposure}, which
#'   will produce population time plots
#' @seealso \code{\link{plot.popTime}}, \code{\link{plot.popTimeExposure}}
#'
#' @import data.table
#' @export


popTime <- function(data, time, event, censored.indicator = NULL,
                    exposure = NULL){

    #data <- DTsim;censored.indicator = NULL;exposure = "z" ; time = "time"; event = "event"
    #names(data)
    varNames <- checkArgsTimeEvent(data = data, time = time, event = event)
    #varNames <- checkArgsTimeEvent(data)
    #varNames <- checkArgsTimeEvent(data, time = time)

    DT <- data.table::as.data.table(data)

    if (is.null(exposure)) {

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
        nLevels <- modifiedEvent$nLevels

        # people with
        # short values of t at the top
        DT[ DT[,order(time)], ycoord:=(nobs:1)]

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
                                        function(i) DT[time >=i & event != 1, .N])]

        # if the 50th percentile number of available subjects at any given
        # point is less than 10, then sample regardless of case status
        if (DT[,quantile(n_available, probs = 0.5)] < 10) {
            message("Sampling from all remaining individuals under study,
                    regardless of event status")
            DT[event == 1,
               n_available := sapply(time,
                                     function(i) DT[time >=i , .N ])]

            DT[event == 1 & n_available > 0,
               yc := sapply(time,
                            function(i)
                                sample(DT[time >= i, ycoord], 1))]

            # use original coordinate if there is no one left to sample from
            DT[event == 1 & n_available == 0, yc := ycoord]
        } else {

            message("Sampling only from individuals who never experienced
                    the event of interest")

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
                   `event status`= checkArgsEventIndicator(data = i, event = varNames$event,
                              censored.indicator = censored.indicator)$event.factor
                        )
                    }
        )

        lapply(l, function(i) {
            nobs <- nrow(i)
            i[ i[,order(time)], ycoord:=(nobs:1)]
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
                                                 function(i) K[time >=i & event != 1, .N])]

            # if the 50th percentile number of available subjects at any given
            # point is less than 10, then sample regardless of case status
            if (K[,quantile(n_available, probs = 0.5)] < 10) {
                message("Sampling from all remaining individuals under study,
                    regardless of event status")
                K[event == 1,
                   n_available := sapply(time,
                                         function(i) K[time >=i , .N ])]

                K[event == 1 & n_available > 0,
                   yc := sapply(time,
                                function(i)
                                    sample(K[time >= i, ycoord], 1))]

                # use original coordinate if there is no one left to sample from
                K[event == 1 & n_available == 0, yc := ycoord]
            } else {

                message("Sampling only from individuals who never experienced
                    the event of interest")

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

        return(lkj)

    }

}





#' @rdname popTime
#' @export
checkArgsTimeEvent <- function(data, time, event) {

    if (missing(time)) {
        if (any(grepl("[\\s\\W_]+time|^time\\b", names(data),
                      ignore.case = TRUE, perl = TRUE))) {
            time <- grep("[\\s\\W_]+time|^time\\b", names(data),
                         ignore.case = TRUE, value = TRUE, perl = TRUE)
            if (length(time)>1)
                warning(paste0("The following variables for time were found in
                              the data: ",paste0(time, collapse = ", "),". '", time[1],
                              "' will be used as the time variable" )) else
                                  print(paste0("'",time,"'",
                                              " will be used as the time variable"),
                                        quote = FALSE)
        } else {
            stop("data does not contain time variable")
        }
    }

    if (missing(event)) {
        if (any(grepl("[\\s\\W_]+event|^event\\b|[\\s\\W_]+status|^status\\b",
                      names(data)[-which(colnames(data)==time[1])],
                      ignore.case = TRUE, perl = TRUE))) {
            event <- grep("[\\s\\W_]+event|^event\\b|[\\s\\W_]+status|^status\\b",
                          names(data)[-which(colnames(data)==time[1])],
                          ignore.case = TRUE, value = TRUE, perl = TRUE)
            if (length(event)>1)
                warning(paste0("The following variables for event were found in
                              the data: ",paste0(event, collapse = ", "),". '", event[1],
                              "' will be used as the event variable" )) else
                                  print(paste0("'",event,"'",
                                              " will be used as the event variable"),
                                        quote = FALSE)
        } else {
            stop("data does not contain event or status variable")
        }
    }

    if (!all(c(time, event) %in% colnames(data))) {
        stop("data does not contain supplied time and/or event variables")
    }

    return(list(time = time[1], event = event[1]))
}




#' Check that Event is in Correct Format
#'
#' Checks for event categories and gives a warning message indicating which
#' level is assumed to be the reference level.
#'
#' @inheritParams popTime
#' @return A list of length two. The first element is the factored event, and
#'   the second element is the numeric representation of the event
#'
#' @export
#' @examples
#'
#' \dontrun{
#' library(survival) # for veteran data
#' checkArgsEventIndicator(data = veteran, event = "celltype", censored.indicator = "smallcell")
#' checkArgsEventIndicator(data = veteran, event = "status")
#' checkArgsEventIndicator(data = veteran, event = "trt") # returns error
#'
#' bmt <- read.csv("https://raw.githubusercontent.com/sahirbhatnagar/casebase/master/inst/extdata/bmtcrr.csv")
#' checkArgsEventIndicator(data = bmt, event = "Sex", censored.indicator = "M")
#' checkArgsEventIndicator(data = bmt, event = "D", censored.indicator = "AML")
#' checkArgsEventIndicator(data = bmt, event = "D", censored.indicator = "AMLL") #returns error
#' checkArgsEventIndicator(data = bmt, event = "Source")
#' checkArgsEventIndicator(data = bmt, event = "Status")
#' checkArgsEventIndicator(data = bmt, event = "Status", censored.indicator = 3)
#' }
#'
checkArgsEventIndicator <- function(data, event, censored.indicator = NULL) {

    isFactor <- is.factor(data[,event])
    isNumeric <- is.numeric(data[,event])
    isCharacter <- is.character(data[,event])

    if (!any(isFactor, isNumeric, isCharacter))
        stop(strwrap("event variable must be either a factor,
                     numeric or character variable", width = 60))

    nLevels <- nlevels(factor(data[,event]))
    if (nLevels < 2) stop(strwrap("event variable must have at least two unique values"))

    if (is.null(censored.indicator)) {

        if (isFactor) {
            slev <- levels(data[,event])
            warning(paste0("censor.indicator not specified. assuming ",
                           slev[1], " represents a censored observation and ",
                           slev[2], " is the event of interest"))
            event.factored <- data[,event]
        }

        if (isCharacter) {
            event.factored <- factor(data[,event])
            slev <- levels(event.factored)
            warning(paste0("censor.indicator not specified. assuming ",
                           slev[1], " represents a censored observation and ",
                           slev[2], " is the event of interest"))
        }

        if (isNumeric) {

            slev <- sort(unique(data[,event]))
            if (!any(slev %in% 0)) stop(strwrap("event is a numeric variable that
                                        doesn't contain 0. if event is a numeric
                                        it must contain some 0's
                                        to indicate censored observations"))
            event.factored <- if (nLevels == 2) factor(data[,event],
                                     labels = c("censored","event")) else
                                         factor(data[,event],
                                                labels = c("censored","event",
                                        paste0("competing event",
                                              if (nLevels >= 4) 1:(nLevels-2))))
        }

    } else {

        if (!(censored.indicator %in% data[,event]) & any(isCharacter, isFactor))
            stop(strwrap("censored.indicator not found in event variable of data"))

        if (isNumeric) {
            warning(strwrap("censored.indicator specified but ignored because
                                event is a numeric variable"))
            slev <- sort(unique(data[,event]))
            if (!any(slev %in% 0)) stop(strwrap("event is a numeric variable that
                                        doesn't contain 0. if event is a numeric
                                        it must contain some 0's
                                        to indicate censored observations"))
            event.factored <- if (nLevels == 2) factor(data[,event],
                                     labels = c("censored","event")) else
                                        factor(data[,event],
                                              labels = c("censored","event",
                                                         paste0("competing event",
                                                if (nLevels >= 4) 1:(nLevels-2))))

        }

        if (isFactor | isCharacter) {

            event.factored <- relevel(factor(data[,event]), censored.indicator)
            slev <- levels(event.factored)
            message(paste0("assuming ",
                           slev[1], " represents a censored observation and ",
                           slev[2], " is the event of interest"))
        }
    }

    return(list(event.factored = event.factored,
                event.numeric = as.numeric((event.factored))-1,
                nLevels = nLevels))

}









# not used ----------------------------------------------------------------

rm(list=ls())


# put the e=0 people at the bottom of the population-time plot and people with
# short values of t at the top
DTsim[ DTsim[,order(time)], ycoord:=(nobs:1)]

# sample y coordinates for each event, so that we can see the incidence density
# on population-time plots. Sampling from people with e=0 or 2 who have an
# observed time t greater than that of a fixed individual who had the event
DT[, yc := 0L]
DT[event == 1, yc := sapply(time,
                            function(i)
                                sample(DT[time >= i & event != 1,ycoord],1) )]

op <- par(mar=c(4.5,4.5,1,1), lwd=2)
with(DT,plot(1,1,type='n', xlim=c(0,end_of_study_time[1]), ylim=c(1,nobs), xlab='Follow-up years', ylab='Population'))
with(DT,segments(rep(0.0, nobs), ycoord, time, ycoord, col='gray'))
with(DT[event==1],points(time, yc, pch=20, col='red', cex=0.5))
par(op)


ggplot(DT, aes(x=0, xend=time, y=ycoord, yend=ycoord)) +
    xlab("Follow-up years") +
    ylab("Population") +
    theme_classic() +
    geom_segment(size=3, colour="grey90") +
    geom_point(aes(x=time, y=yc), data = DT[event==1], size=1, colour="red") +
    theme(axis.text=element_text(size=12, face='bold')) +
    labs(title='Population time plot. Red dots are cases.')




rm(list=ls())
library(survival)
data(veteran)
veteran$prior <- factor(veteran$prior, levels = c(0, 10))
veteran$celltype <- factor(veteran$celltype,
                           levels = c('large', 'squamous', 'smallcell', 'adeno'))
veteran$trt <- factor(veteran$trt, levels = c(1, 2))
isEmpty <- function(x) {
    return(length(x) == 0)
}



head(veteran)
rm(DTv)
DTv <- as.data.table(veteran)
nobs <- nrow(DTv)

DTv[,table(status)]

DTv[ DTv[,order(time)], ycoord:=(nobs:1)]

# sample y coordinates for each event, so that we can see the incidence density
# on population-time plots. Sampling from people with e=0 or 2 who have an
# observed time t greater than that of a fixed individual who had the event
DTv[, yc := 0L]
str(DTv)

# need to find out how many controls are left to sample from
DTv[, n_available := 0L]
# DTv[status == 1, n_available := sapply(time, function(i) DTv[time >=i & status != 1, .N ] )]

DTv[status == 1, n_available := sapply(time, function(i) DTv[time >=i, .N ] )]

# DTv[status == 1 & n_available > 0,
#     yc := sapply(time,
#                  function(i)
#                      sample(DTv[time >= i & status != 1,ycoord],1) )]

DTv[status == 1 & n_available > 0,
    yc := sapply(time,
                 function(i)
                     sample(DTv[time >= i,ycoord],1) )]


ggplot(DTv, aes(x=0, xend=time, y=ycoord, yend=ycoord)) +
    xlab("Follow-up years") +
    ylab("Population") +
    theme_classic() +
    geom_segment(size=3, colour="grey90") +
    geom_point(aes(x=time, y=yc), data = DTv[status==1], size=1, colour="red") +
    theme(axis.text=element_text(size=12, face='bold')) +
    labs(title='Population time plot. Red dots are cases.')


with(DTv, Surv(time, status))




str(veteran)

Surv()





DT <- read.csv("https://raw.githubusercontent.com/sahirbhatnagar/casebase/master/inst/extdata/bmtcrr.csv")
nobs <- nrow(DT)
ftime <- DT$ftime
ord <- order(ftime, decreasing=TRUE)
plot(0, type='n', xlim=c(0, max(ftime)), ylim=c(0, nobs),
     xlab='Follow-up time', ylab='Population')
#segments(rep(0.0, nobs), 1:nobs, ftime[ord], 1:nobs, col='gray25')
polygon(c(0, max(ftime), ftime[ord], 0), c(0, 0, 1:nobs, nobs), col = "gray90")
cases <- DT$Status %in% c(1, 2)
colour <- c("red", "blue")[DT$Status[cases]]









head(veteran)
names(veteran)
dev.off()
data = as.data.table(veteran)
str(data)

ggplot(veteran, aes(time = time, event = status)) +
    stat_poptime(na.rm = FALSE)

veteran$status






# get row number of subjects who have an event==1, and either covariate value
cond <- DT[e==1 & z %in% c(0,1), which=T]

# get event time of subjects who have an event==1, and either covariate value
etimes <- DT[e==1 & z %in% c(0,1), t]

