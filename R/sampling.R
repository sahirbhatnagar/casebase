# This is where all functions related to sampling should appear
expit <- function(x) 1/(1 + exp(-x))
logit <- function(p) log(p)-log(1-p)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

# R function to create case-base dataset for use in fitting
# parametric hazard functions via logistic regression
# see Hanley and Miettinen, Int J Biostatistics 2009

#' @param ds source dataset
#' @param event.var event variable (1=event)
#' @param t.var event (or censoring) time
#' @param i.var intervention (tx) variable
#' @param id.var patient identifier
#' @param x.vars vector of names of regressor variables
#' @param b.c.ratio (integer) ratio, size of base series : case series
#' @param random if TRUE,
#' @return returns dataset with b+c rows of person-moments (p.m), x.vars,
#'   offset, and an indicator variable y (1 if p.m represents an event, 0
#'   otherwise)

sampleCaseBase <- function(ds, event.var, t.var,
                           i.var, id.var, x.vars,
                           b.c.ratio, random = TRUE) {

    n.subjects <- length(ds[,t.var]) # no. of subjects
    B <- sum(ds[,t.var])             # total person-time in base
    c <- sum(ds[,event.var])          # no. of cases (events)
    b <- b.c.ratio * c               # size of base series
    offset <- log(B / b)            # offset so intercept = log(ID | x, t = 0 )

    if (random) {
        p <- ds[,t.var]/B
        who <- sample(n.subjects, b, replace = TRUE, prob = p)
        b.series <- ds[who,]
        b.series <- b.series[,c(i.var,id.var,x.vars,t.var)]
        b.series$y <- 0
        b.series[,t.var] <- runif(b)*b.series[,t.var]
        b.series$o <- offset
    }

    if (!random) {
        d.t <- B/(b+1)
        p.sum <- c(0) #Allocate memory first!!
        for (i in 1:n.subjects) {
            p.sum <- c(p.sum, p.sum[i] + ds[i,t.var])
        }
        every.d.t <- B*(1:b)/(b+1)
        who <- findInterval(every.d.t, p.sum)
        #print(who)
        b.series <- ds[who,]
        b.series <- b.series[, c(i.var,id.var,x.vars, t.var)]
        b.series$y <- 0
        b.series[,t.var] <- every.d.t - p.sum[who]
        b.series$o <- offset
    }

    c.series <- ds[ ds[event.var]==1, ]
    c.series <- c.series[,c(i.var, id.var, x.vars, t.var)]
    c.series$y <- 1
    c.series[, t.var] <- c.series[, t.var]
    c.series$o <- offset
    c.b.series <- rbind(c.series, b.series)

    return(c.b.series)
}
