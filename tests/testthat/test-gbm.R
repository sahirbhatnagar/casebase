context("GBMs")

# Skip tests if gbm is not installed
testthat::skip_if_not_installed("gbm")

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

formula_gbm <- formula(paste(c("event ~ ftime", "Z",
                               paste0("V", 1:10)),
                             collapse = " + "))

# Fitting----
test_that("no error in fitting gbm", {
    fitDF <- try(fitSmoothHazard(formula_gbm, data = DF_ext, time = "ftime",
                                 family = "gbm"),
                 silent = TRUE)
    fitDT <- try(fitSmoothHazard(formula_gbm, data = DT_ext, time = "ftime",
                                 family = "gbm"),
                 silent = TRUE)

    expect_false(inherits(fitDF, "try-error"))
    expect_false(inherits(fitDT, "try-error"))
})

test_that("warnings witn non-linear functions of time or interactions", {
    expect_warning(fitSmoothHazard(event ~ log(ftime) + Z,
                                   data = DF, time = "ftime", family = "gbm"),
                   regexp = "gbm may throw an error")
    expect_warning(fitSmoothHazard(event ~ ftime * Z,
                                   data = DF, time = "ftime", family = "gbm"),
                   regexp = "gbm may throw an error")
})

# Absolute risk----
fitDF_gbm <- fitSmoothHazard(event ~ ftime + Z, data = DF, time = "ftime",
                             family = "gbm", ratio = 10)
fitDT_gbm <- fitSmoothHazard(event ~ ftime + Z, data = DT, time = "ftime",
                             family = "gbm", ratio = 10)

newDT <- data.table("Z" = c(0, 1))
newDF <- data.frame("Z" = c(0, 1))

test_that("no error in fitting gbm", {
    riskDF <- try(absoluteRisk(fitDF_gbm, time = 0.5, newdata = newDF,
                               n.trees = 100, nsamp = 500),
                  silent = TRUE)
    riskDT <- try(absoluteRisk(fitDT_gbm, time = 0.5, newdata = newDT,
                               n.trees = 100, nsamp = 500),
                  silent = TRUE)
    riskDF_mc <- try(absoluteRisk(fitDF_gbm, time = 0.5, newdata = newDF,
                                  n.trees = 100, nsamp = 10,
                                  method = "montecarlo"),
                     silent = TRUE)
    riskDT_mc <- try(absoluteRisk(fitDT_gbm, time = 0.5, newdata = newDT,
                                  n.trees = 100, nsamp = 10,
                                  method = "montecarlo"),
                     silent = TRUE)

    expect_false(inherits(riskDF, "try-error"))
    expect_false(inherits(riskDT, "try-error"))
    expect_false(inherits(riskDF_mc, "try-error"))
    expect_false(inherits(riskDT_mc, "try-error"))
})

test_that("should compute risk when time and newdata aren't provided", {
    # To pass the test, I had to increase nsamp
    # This may have something to do with the way we use gbm
    # or it could be that the estimated hazard is highly non-smooth
    # In any case, we will have to test gbm more to see what's going on.
    skip_on_cran()
    fitDF_gbm_red <- fitDF_gbm
    fitDF_gbm_red$originalData <- fitDF_gbm$originalData[c(1:5, 101:105), ]
    absRiskDF_gbm <- absoluteRisk(fitDF_gbm_red, n.trees = 100, nsamp = 500)

    fitDT_gbm_red <- fitDT_gbm
    fitDT_gbm_red$originalData <- fitDT_gbm$originalData[c(1:5, 101:105), ]
    absRiskDT_gbm <- absoluteRisk(fitDT_gbm_red, n.trees = 100, nsamp = 500)

    expect_true("risk" %in% names(absRiskDF_gbm))
    expect_true("risk" %in% names(absRiskDT_gbm))
})

test_that("output probabilities", {
    riskDF_gbm <- absoluteRisk(fitDF_gbm, time = 0.5, newdata = newDF,
                               family = "gbm", n.trees = 100, nsamp = 500)
    riskDT_gbm <- absoluteRisk(fitDT_gbm, time = 0.5, newdata = newDT,
                               family = "gbm", n.trees = 100, nsamp = 500)

    expect_true(all(riskDF_gbm >= 0))
    expect_true(all(riskDT_gbm >= 0))
    expect_true(all(riskDF_gbm <= 1))
    expect_true(all(riskDT_gbm <= 1))
})

# Summary method
test_that("no error in summary method for gbm", {
    sumDF <- try(print(summary(fitDF_gbm)),
                 silent = TRUE)
    sumDT <- try(print(summary(fitDT_gbm)),
                 silent = TRUE)

    expect_false(inherits(sumDF, "try-error"))
    expect_false(inherits(sumDT, "try-error"))
})

# Matrix interface----
N <- 1000; p <- 30
nzc <- 0.33 * p
x <- matrix(rnorm(N * p), N, p)
dimnames(x)[[2]] <- paste0("x", seq_len(p))
beta <- rnorm(nzc)
fx <- x[, seq(nzc)] %*% (0.33 * beta)
hx <- exp(fx)
ty <- rexp(N, hx)
tcens <- rbinom(n = N,
                prob = 0.3,
                size = 1) # censoring indicator
y <- cbind(time = ty, status = 1 - tcens) # y=Surv(ty,1-tcens) with survival

muffler <- function(msg) {
    if (any(grepl("condition has length > 1", msg))) {
        invokeRestart("muffleWarning")
        }
}

skip_next_tests <- (Sys.getenv("_R_CHECK_LENGTH_1_CONDITION_") == "true" ||
                        Sys.getenv("_R_CHECK_LENGTH_1_LOGIC2_") == "true")

testthat::skip_if(skip_next_tests,
                  "gbm throws an error because it checks for equality of class\ninstead of using inherits (version 2.1.8)")

test_that("no error in fitting fitSmoothHazard.fit", {
    # gbm throws a warning because it checks for equality of class
    # instead of using inherits (version 2.1.8)
    fit_gbm <- try(withCallingHandlers(fitSmoothHazard.fit(x, y, time = "time",
                                                           event = "status",
                                       family = "gbm", ratio = 10),
                                       warning = muffler),
                   silent = TRUE)

    expect_false(inherits(fit_gbm, "try-error"))
})

test_that("warnings witn non-linear functions of time or interactions", {
    skip_if_not_installed("splines")
    library(splines)
    expect_warning(fitSmoothHazard.fit(x, y, formula_time = ~ log(time),
                                       time = "time", event = "status",
                                       family = "gbm", ratio = 10),
                   regexp = "gbm may throw an error")
    expect_warning(fitSmoothHazard.fit(x, y, formula_time = ~ bs(time),
                                       time = "time", event = "status",
                                       family = "gbm", ratio = 10),
                   regexp = "gbm may throw an error")
})
