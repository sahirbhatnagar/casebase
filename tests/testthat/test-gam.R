context("GAMs")

# Skip tests if mgcv is not installed
testthat::skip_if_not_installed("mgcv")

# Create data----
n <- 100
alp <- 0.05
lambda_t0 <- 1
lambda_t1 <- 3

times <- c(rexp(n = n, rate = lambda_t0),
           rexp(n = n, rate = lambda_t1))
censor <- rexp(n = 2 * n, rate = -log(alp))

times_c <- pmin(times, censor)
event_c <- 1 * (times < censor)

DF <- data.frame("ftime" = times_c,
                 "event" = event_c,
                 "Z" = c(rep(0, n),
                         rep(1, n)))
DT <- data.table("ftime" = times_c,
                 "event" = event_c,
                 "Z" = c(rep(0, n),
                         rep(1, n)))

extra_vars <- matrix(rnorm(10 * n), ncol = 10)
DF_ext <- cbind(DF, as.data.frame(extra_vars))
DT_ext <- cbind(DT, as.data.table(extra_vars))

formula_gam <- formula(paste(c("event ~ s(ftime)", "Z",
                               paste0("V", 1:10)),
                             collapse = " + "))

# Fitting----
test_that("no error in fitting gam", {
    fitDF <- try(fitSmoothHazard(formula_gam, data = DF_ext, time = "ftime",
                                 family = "gam"),
                 silent = TRUE)
    fitDT <- try(fitSmoothHazard(formula_gam, data = DT_ext, time = "ftime",
                                 family = "gam"),
                 silent = TRUE)

    expect_false(inherits(fitDF, "try-error"))
    expect_false(inherits(fitDT, "try-error"))
})

# Absolute risk----
fitDF_gam <- fitSmoothHazard(event ~ s(ftime) + Z, data = DF,
                             time = "ftime", family = "gam", ratio = 10)
fitDT_gam <- fitSmoothHazard(event ~ s(ftime) + Z, data = DT,
                             time = "ftime", family = "gam", ratio = 10)

newDT <- data.table("Z" = c(0, 1))
newDF <- data.frame("Z" = c(0, 1))

test_that("no error in abs risk for gam", {
    riskDF <- try(absoluteRisk(fitDF_gam, time = 0.5, newdata = newDF),
                  silent = TRUE)
    riskDT <- try(absoluteRisk(fitDT_gam, time = 0.5, newdata = newDT),
                  silent = TRUE)
    riskDF_mc <- try(absoluteRisk(fitDF_gam, time = 0.5, newdata = newDF,
                                  method = "montecarlo"),
                     silent = TRUE)
    riskDT_mc <- try(absoluteRisk(fitDT_gam, time = 0.5, newdata = newDT,
                                  method = "montecarlo"),
                     silent = TRUE)

    expect_false(inherits(riskDF, "try-error"))
    expect_false(inherits(riskDT, "try-error"))
    expect_false(inherits(riskDF_mc, "try-error"))
    expect_false(inherits(riskDT_mc, "try-error"))
})

test_that("should compute risk when time and newdata aren't provided", {
    absRiskDF_gam <- absoluteRisk(fitDF_gam)
    absRiskDT_gam <- absoluteRisk(fitDT_gam)

    expect_true("risk" %in% names(absRiskDF_gam))
    expect_true("risk" %in% names(absRiskDT_gam))
})

test_that("output probabilities", {
    riskDF_gam <- absoluteRisk(fitDF_gam, time = 0.5, newdata = newDF,
                               family = "gam")
    riskDT_gam <- absoluteRisk(fitDT_gam, time = 0.5, newdata = newDT,
                               family = "gam")

    expect_true(all(riskDF_gam >= 0))
    expect_true(all(riskDT_gam >= 0))
    expect_true(all(riskDF_gam <= 1))
    expect_true(all(riskDT_gam <= 1))
})

# Summary method
test_that("no error in summary method for gam", {
    sumDF <- try(print(summary(fitDF_gam)),
                  silent = TRUE)
    sumDT <- try(print(summary(fitDT_gam)),
                  silent = TRUE)

    expect_false(inherits(sumDF, "try-error"))
    expect_false(inherits(sumDT, "try-error"))
})
