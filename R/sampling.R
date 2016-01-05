#' Create case-base dataset for use in fitting parametric hazard functions
#'
#' This function implements the case-base sampling approach described in Hanley
#' and Miettinen, Int J Biostatistics 2009. It can be used to fit smooth-in-time
#' parametric functions easily, via logistic regression.
#'
#' It is assumed that \code{data} contains two columns named \code{time} and
#' \code{event} (see the function \code{\link{Surv}} in the package
#' \code{survival}).
#'
#' @param data The source dataset; see Details.
#' @param ratio Integer, giving the ratio of the size of the base series to that
#'   of the case series. Defaults to 10.
#' @param type There are currently two sampling procedures available:
#'   \code{uniform}, where person-moments are sampled uniformly across
#'   individuals and follow-up time; and \code{multinomial}, where individuals
#'   are sampled proportionally to their follow-up time.
#' @return The function returns a dataset, with the same format as the source
#'   dataset, and where each row corresponds to a person-moment sampled from the
#'   case or the base series. otherwise)
sampleCaseBase <- function(data, ratio = 10, type = c("uniform", "multinomial")) {

    type <- match.arg(type)
    # Create survival object from dataset
    if (! all(c("time", "event") %in% colnames(data))) {
        stop("data should contain two columns named time and event, respectively.")
    }
    survObj <- with(data, survival::Surv(time, event))

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
        bSeries <- as.matrix(survObj[who, ])
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
        pSum <- c(0, cumsum(survObj[, "time"]))
        everyDt <- B*(1:b)/(b+1)
        who <- findInterval(everyDt, pSum)
        bSeries <- as.matrix(survObj[who, ])
        bSeries[, "status"] <- 0
        bSeries[, "time"] <- everyDt - pSum[who]
    }

    # Next commented line will break on data.table
    # bSeries <- cbind(bSeries, data[who, colnames(data) != c("time", "event")])
    bSeries <- cbind(bSeries, subset(data, select = (colnames(data) != c("time", "event")))[who,])
    names(bSeries)[names(bSeries) == "status"] <- "event"

    cSeries <- data[data$event == 1,]
    # cSeries <- survObj[survObj[, "status"] == 1, ]
    # cSeries <- cSeries[, c(i.var, id.var, x.vars, time)]
    # cSeries$y <- 1
    # cSeries[, time] <- cSeries[, time]

    # Combine case and base series
    cbSeries <- rbind(cSeries, bSeries)
    cbSeries <- data.table::data.table(cbSeries, offset = rep_len(offset, nrow(cbSeries)))

    class(cbSeries) <- c("cbData", class(cbSeries))
    return(cbSeries)
}
