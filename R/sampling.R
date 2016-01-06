#' Create case-base dataset for use in fitting parametric hazard functions
#'
#' This function implements the case-base sampling approach described in Hanley
#' and Miettinen, Int J Biostatistics 2009. It can be used to fit smooth-in-time
#' parametric functions easily, via logistic regression.
#'
#' It is assumed that \code{data} contains the two columns corresponding to the
#' supplied time and event variables. If either the \code{time} or \code{event}
#' argument is missing, the function looks for columns named \code{"time"},
#' \code{"event"}, or \code{"status"}.
#'
#' @param data a data.frame or data.table containing the source dataset.
#' @param time a character string giving the name of the time variable. See
#'   Details.
#' @param event a character string giving the name of the event variable. See
#'   Details.
#' @param ratio Integer, giving the ratio of the size of the base series to that
#'   of the case series. Defaults to 10.
#' @param type There are currently two sampling procedures available:
#'   \code{uniform}, where person-moments are sampled uniformly across
#'   individuals and follow-up time; and \code{multinomial}, where individuals
#'   are sampled proportionally to their follow-up time.
#' @return The function returns a dataset, with the same format as the source
#'   dataset, and where each row corresponds to a person-moment sampled from the
#'   case or the base series. otherwise)
sampleCaseBase <- function(data, time, event, ratio = 10, type = c("uniform", "multinomial")) {
    if (missing(time)) {
        if ("time" %in% colnames(data)) {
            time <- "time"
        } else {
            stop("data does not contain time variable")
        }
    }
    if (missing(event)) {
        if ("event" %in% colnames(data)) {
            event <- "event"
        } else {
            if ("status" %in% colnames(data)) {
                event <- "status"
            } else {
                stop("data does not contain event or status variable")
            }
        }
    }
    if (!all(c(time, event) %in% colnames(data))) {
        stop("data does not contain supplied time and event variables")
    }
    type <- match.arg(type)
    # Create survival object from dataset
    survObj <- survival::Surv(subset(data, select=(names(data) == time)),
                              subset(data, select=(names(data) == event)))

    n <- nrow(survObj) # no. of subjects
    B <- sum(survObj[, "time"])             # total person-time in base
    c <- sum(survObj[, "status"])          # no. of cases (events)
    b <- ratio * c               # size of base series
    offset <- log(B / b)            # offset so intercept = log(ID | x, t = 0 )

    if (type == "uniform") {
        # The idea here is to sample b individuals, with replacement, and then
        # to sample uniformly a time point for each of them. The sampled time
        # point must lie between the beginning and the end of follow-up
        p <- survObj[, "time"]/B
        who <- sample(n, b, replace = TRUE, prob = p)
        bSeries <- survObj[who, ]
        bSeries[, "status"] <- 0
        bSeries[, "time"] <- runif(b) * bSeries[, "time"]
    }

    if (type == "multinomial") {
        # Multinomial sampling: probability of individual contributing a
        # person-moment to base series is proportional to time variable
        # dt <- B/(b+1)
        # pSum <- c(0) #Allocate memory first!!
        # for (i in 1:n) {
        #     pSum <- c(pSum, pSum[i] + survObj[i, "time"])
        # }
        pSum <- c(0, cumSum(survObj[, "time"]))
        everyDt <- B*(1:b)/(b+1)
        who <- findInterval(everyDt, pSum)
        bSeries <- survObj[who, ]
        bSeries[, "status"] <- 0
        bSeries[, "time"] <- everyDt - pSum[who]
    }

    # Next commented line will break on data.table
    # bSeries <- cbind(bSeries, data[who, colnames(data) != c("time", "event")])
    bSeries <- cbind(bSeries, subset(data, select = (colnames(data) != c(time, event)))[who,])
    names(bSeries)[names(bSeries) == "status"] <- event

    cSeries <- data[which(subset(data, select=(names(data) == event)) == 1),]
    # cSeries <- survObj[survObj[, "status"] == 1, ]
    # cSeries <- cSeries[, c(i.var, id.var, x.vars, time)]
    # cSeries$y <- 1
    # cSeries[, time] <- cSeries[, time]

    # Combine case and base series
    cbSeries <- rbind(cSeries, bSeries)
    cbSeries <- cbind(cbSeries, rep_len(offset, nrow(cbSeries)))
    names(cbSeries)[ncol(cbSeries)] <- "offset"

    class(cbSeries) <- c("cbData", class(cbSeries))
    return(cbSeries)
}
