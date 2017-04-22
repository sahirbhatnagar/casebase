context("popTime methods")

nobs <- 500

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
DTsim <- data.table(ID = seq_len(nobs), z = rbinom(nobs, 1, 0.5))
setkey(DTsim, ID)
DTsim[,`:=`(event_time = rweibull(nobs, a1, b1 * exp(z * c1)^(-1/a1)),
             competing_time = rweibull(nobs, a2, b2 * exp(z * c2)^(-1/a2)),
             end_of_study_time = eost)]
DTsim[,`:=`(event = 1 * (event_time < competing_time) +
                2 * (event_time >= competing_time),
            time = pmin(event_time, competing_time))]
DTsim[time >= end_of_study_time, event := 0]
DTsim[time >= end_of_study_time, time := end_of_study_time]

DFsim <- data.frame("z" = DTsim[,z],
                    "event" = DTsim[,event],
                    "time" = DTsim[,time])

test_that("no error in popTime with data.frame or data.table", {
    out1 <- try(popTime(data = DTsim, time = "time", event = "event"))
    out2 <- try(popTime(data = DFsim, time = "time", event = "event"))

    expect_false(inherits(out1, "try-error"))
    expect_false(inherits(out2, "try-error"))
})

test_that("no error in stratified popTime with data.frame or data.table", {
    out1 <- try(popTime(data = DTsim, time = "time",
                        event = "event", exposure = "z"))
    out2 <- try(popTime(data = DFsim, time = "time",
                        event = "event", exposure = "z"))

    expect_false(inherits(out1, "try-error"))
    expect_false(inherits(out2, "try-error"))
})
