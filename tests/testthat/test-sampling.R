context("Sampling")

test_that("no error in sampling with data.frame or data.table", {
    nobs <- 5000
    tlim <- 10

    # simulation parameters
    b1 <- 200
    b2 <- 50

    # event type 0-censored, 1-event of interest, 2-competing event
    # t observed time/endpoint
    # z is a binary covariate
    DT <- data.table(z = rbinom(nobs, 1, 0.5))
    DT[,`:=`("t_event" = rweibull(nobs, 1, b1),
              "t_comp" = rweibull(nobs, 1, b2))]
    DT[,`:=`("event" = 1 * (t_event < t_comp) + 2 * (t_event >= t_comp),
             "time" = pmin(t_event, t_comp))]
    DT[time >= tlim, `:=`("event" = 0, "time" = tlim)]

    DF <- data.frame(z = rbinom(nobs, 1, 0.5),
                     t_event = rweibull(nobs, 1, b1),
                     t_comp = rweibull(nobs, 1, b2))
    DF$event <- with(DF, 1 * (t_event < t_comp) + 2 * (t_event >= t_comp))
    DF$time <- with(DF, pmin(t_event, t_comp))
    DF[DF$time >= tlim, ]$event <- 0
    DF[DF$time >= tlim, ]$time <- tlim

    out1 <- try(sampleCaseBase(DT, time = "time", event = "event", comprisk = TRUE))
    out2 <- try(sampleCaseBase(DF, time = "time", event = "event", comprisk = TRUE))

    expect_false(inherits(out1, "try-error"))
    expect_false(inherits(out2, "try-error"))
})
