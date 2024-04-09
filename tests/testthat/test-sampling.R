context("Sampling")
set.seed(12345)

# CRAN skip atlas check fix
testthat::skip_if(grepl(pattern = "atlas", sessionInfo()$BLAS,
                        ignore.case = TRUE))

# Create simulated data with competing risks----
nobs <- 500
tlim <- 10
b1 <- 200
b2 <- 50

# event type 0-censored, 1-event of interest, 2-competing event
# t observed time/endpoint
# z is a binary covariate
DT <- data.table(z = rbinom(nobs, 1, 0.5))
DT[, `:=`("t_event" = rweibull(nobs, 1, b1),
         "t_comp" = rweibull(nobs, 1, b2))]
DT[, `:=`("event" = 1 * (t_event < t_comp) + 2 * (t_event >= t_comp),
         "time" = pmin(t_event, t_comp))]
DT[time >= tlim, `:=`("event" = 0, "time" = tlim)]
DT[, c("t_event", "t_comp") := NULL]

DF <- data.frame(z = rbinom(nobs, 1, 0.5),
                 t_event = rweibull(nobs, 1, b1),
                 t_comp = rweibull(nobs, 1, b2))
DF$event <- with(DF, 1 * (t_event < t_comp) + 2 * (t_event >= t_comp))
DF$time <- with(DF, pmin(t_event, t_comp))
DF[DF$time >= tlim, ]$event <- 0
DF[DF$time >= tlim, ]$time <- tlim
DF$t_event <- NULL
DF$t_comp <- NULL

test_that("Expect error with competing risk but compRisk is not specified", {
    expect_error(sampleCaseBase(DT, time = "time", event = "event"))
    expect_error(sampleCaseBase(DF, time = "time", event = "event"))
})

test_that("no error in sampling with data.frame or data.table", {
    out1 <- try(sampleCaseBase(DT, time = "time", event = "event",
                               comprisk = TRUE))
    out2 <- try(sampleCaseBase(DF, time = "time", event = "event",
                               comprisk = TRUE))

    expect_false(inherits(out1, "try-error"))
    expect_false(inherits(out2, "try-error"))
})

test_that("detect time variable within data.frame or data.table", {
    out1 <- try(sampleCaseBase(DT, event = "event", comprisk = TRUE))
    out2 <- try(sampleCaseBase(DF, event = "event", comprisk = TRUE))

    expect_false(inherits(out1, "try-error"))
    expect_false(inherits(out2, "try-error"))
})

test_that("detect event variable within data.frame or data.table", {
    out1 <- try(sampleCaseBase(DT, time = "time", comprisk = TRUE))
    out2 <- try(sampleCaseBase(DF, time = "time", comprisk = TRUE))

    expect_false(inherits(out1, "try-error"))
    expect_false(inherits(out2, "try-error"))
})

# Test checkArgsEventIndicator----
data("bmtcrr") # from casebase
bmtcrr$Sex <- as.character(bmtcrr$Sex)

test_that("no error with different types of event variables", {
    out1 <- try(checkArgsEventIndicator(data = survival::veteran,
                                        event = "celltype",
                                        censored.indicator = "smallcell"))
    out2 <- try(checkArgsEventIndicator(data = survival::veteran,
                                        event = "status"))
    out3 <- try(checkArgsEventIndicator(data = bmtcrr, event = "Sex",
                                        censored.indicator = "M"))

    expect_false(inherits(out1, "try-error"))
    expect_false(inherits(out2, "try-error"))
    expect_false(inherits(out3, "try-error"))
})

test_that("warning when baseline not specified", {
    expect_warning(checkArgsEventIndicator(data = survival::veteran,
                                           event = "celltype"),
                   regexp = "censor.indicator not specified")
    expect_warning(checkArgsEventIndicator(data = bmtcrr, event = "Sex"),
                   regexp = "censor.indicator not specified")
})

test_that("warning when baseline is specified for numeric", {
    expect_warning(checkArgsEventIndicator(data = survival::veteran,
                                           event = "status",
                                           censored.indicator = 0L),
                   regexp = "censored.indicator specified but ignored")
    expect_warning(checkArgsEventIndicator(data = bmtcrr, event = "Status",
                                           censored.indicator = 0L),
                   regexp = "censored.indicator specified but ignored")
})
