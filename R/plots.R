library(data.table)
library(magrittr)
library(ggplot2)

isEmpty <- function(x) {
    return(length(x) == 0)
}

rm(list=ls())
set.seed(1)
nobs <- 5000

# simulation parameters
a1 <- 1.0
b1 <- 200
a2 <- 1.0
b2 <- 50
c1 <- 0.0
c2 <- 0.0

# end of study time
eost <- 10

# e event type 0-censored, 1-event of interest, 2-competing event
# t observed time/endpoint
# z is a binary covariate
DT <- data.table(ID = seq_len(nobs), z=rbinom(nobs, 1, 0.5))

setkey(DT, ID)

DT[,`:=` (event_time = rweibull(nobs, a1, b1 * exp(z * c1)^(-1/a1)),
          competing_time = rweibull(nobs, a2, b2 * exp(z * c2)^(-1/a2)),
          end_of_study_time = eost)]

DT[,`:=`(event = 1 * (event_time < competing_time) + 2 * (event_time >= competing_time),
         time = pmin(event_time, competing_time))]

DT[time >= end_of_study_time, event := 0]

DT[time >= end_of_study_time, time:=end_of_study_time]

# put the e=0 people at the bottom of the population-time plot and people with
# short values of t at the top
DT[ DT[,order(time)], ycoord:=(nobs:1)]

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



#' Population Time Plot
#'
#' Draws a population time plot to give a visual representation of incidence
#' density
#'
#' @param data a data.frame or data.table containing the source dataset.
#' @param time a character string giving the name of the time variable. See
#'   Details.
#' @param event a character string giving the name of the event variable. See
#'   Details.
#' @param title plot title. By default there is no title
#' @param xlab x-axis label. Default is \code{"Follow-up years"}
#' @param ylab y-axis label. Default is \code{"Population"}
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
#'   The following regular expressions are used for the time and event variables:
#'   \describe{
#'   \item{time}{\code{"[\\s\\W_]+time|^time\\b"}}
#'   \item{event}{\code{"[\\s\\W_]+event|^event\\b|[\\s\\W_]+status|^status\\b"}}
#'   }
#'
#'   This allows for \code{"time"} to be preceded or followed by one or more
#'   white space characters, one or more non-word characters or one or more
#'   underscores. For example, the following column names would be recognized by
#'   the function as the \code{"time"} variable: \code{"time of death",
#'   "death_time", "Time", "time", "diagnosis_time", "time.diag", "diag__time"}.
#'   But the following will not be recognized: \code{"diagtime","eventtime",
#'   "Timediag"}
#'
#'
#' @import data.table
#' @export

rm(list = ls())
popTimePlot <- function(data, time, event, title = NULL, xlab = NULL, ylab = NULL){

    ##### Change the code so that you modify the column names of the data.table,
    # then you wont have to deal with variables in Data.table expressions

    # data = veteran;
    # names(veteran)
    # colnames(data)[grep("^time", colnames(data))] <- "time of death"
    # colnames(data)[grep("^status", colnames(data))] <- "death status"
    # colnames(data)[grep("^diagtime", colnames(data))] <- "time status"

    # names(data)
    varNames <- checkArgsTimeEvent(data = data, time = time, event = event)

    DT <- as.data.table(data)
    nobs <- nrow(DT)

    # names(veteran)
    setnames(DT,varNames$time, "time")
    setnames(DT,varNames$event, "event")

    # people with
    # short values of t at the top
    DT[ DT[,order(time)], ycoord:=(nobs:1)]

    # sample y coordinates for each event, so that we can see the incidence density
    # on population-time plots. Sampling from people who have an
    # observed time t greater than that of a fixed individual who had the event
    DT[, yc := 0L]

    # need to find out how many controls are left to sample from
    DT[, n_available := 0L]

    DT[event == 1, n_available := sapply(time, function(i) DT[time >=i, .N ] )]

    DT[event == 1 & n_available > 0,
        yc := sapply(time,
                     function(i)
                         sample(DT[time >= i,ycoord],1) )]

    p <- ggplot(DT, aes(x=0, xend=time, y=ycoord, yend=ycoord))

    p +
        # if(is.null(xlab)) xlab("Follow-up years") else xlab(xlab) +
        # if(is.null(ylab)) ylab("Population") else ylab(ylab) +
        # if(!is.null(title)) labs(title=title) +
        # xlab("Follow-up years") +
        # ylab("Population") +
        theme_classic() +
        geom_segment(size=3, colour="grey90") +
        geom_point(aes(x=time, y=yc), data = DT[event==1], size=1, colour="red") +
        theme(axis.text=element_text(size=12, face='bold'))

}


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
        stop("data does not contain supplied time and event variables")
    }

    return(list(time = time[1], event = event[1]))
}

dev.off()
popTimePlot(veteran)



StatLm <- ggproto("StatLm", Stat,
            required_aes = c("x", "y"),
            compute_group = function(data, scales, params, n = 100, formula = y ~ x) {
                rng <- range(data$x, na.rm = TRUE)
                grid <- data.frame(x = seq(rng[1], rng[2], length = n))

                mod <- lm(formula, data = data)
                grid$y <- predict(mod, newdata = grid)

                grid
            }
)

StatLm <- ggproto("StatLm", Stat,
                  required_aes = c("x", "y"),
                  compute_group = function(data, scales, params, n = 100, formula = y ~ x) {
                      rng <- range(data$x, na.rm = TRUE)
                      grid <- data.frame(x = seq(rng[1], rng[2], length = n))

                      mod <- lm(formula, data = data)
                      grid$y <- predict(mod, newdata = grid)

                      grid
                  }
)


stat_lm <- function(mapping = NULL, data = NULL, geom = "line",
                    position = "identity", na.rm = FALSE, show.legend = NA,
                    inherit.aes = TRUE, n = 50, formula = y ~ x,
                    ...) {
    layer(
        stat = StatLm, data = data, mapping = mapping, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(n = n, formula = formula, na.rm = na.rm, ...)
    )
}




ggplot(mpg, aes(displ, hwy)) +
    geom_point() +
    stat_lm(formula = y ~ poly(x, 10)) +
    stat_lm(formula = y ~ poly(x, 10), geom = "point", colour = "red", n = 20)












# get row number of subjects who have an event==1, and either covariate value
cond <- DT[e==1 & z %in% c(0,1), which=T]

# get event time of subjects who have an event==1, and either covariate value
etimes <- DT[e==1 & z %in% c(0,1), t]

