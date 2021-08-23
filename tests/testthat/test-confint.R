# Handling warning messages coming from montecarlo integration
handler_validmc <- function(msg) {
    if (any(grepl("out of range", msg))) invokeRestart("muffleWarning")
}

# Setup dataset----
n <- 50
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

fitDF <- fitSmoothHazard(event ~ Z, data = DF, time = "ftime", ratio = 10)
fitDT <- fitSmoothHazard(event ~ Z, data = DT, time = "ftime", ratio = 10)

# Start tests----
absDF <- absoluteRisk(fitDF, time = 1, newdata = DF[1, ])
absDT <- absoluteRisk(fitDF, time = 1, newdata = DT[1, ])

test_that("no error in confint with one covariate profile", {
    foo1 <- try(confint(absDF, fitDF, nboot = 10),
                silent = TRUE)
    foo2 <- try(confint(absDT, fitDT, nboot = 10),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})

absDF <- absoluteRisk(fitDF, time = 1, newdata = DF[c(1, n+1), ])
absDT <- absoluteRisk(fitDF, time = 1, newdata = DT[c(1, n+1), ])

test_that("no error in confint with two covariate profiles", {
    foo1 <- try(confint(absDF, fitDF, nboot = 10),
                silent = TRUE)
    foo2 <- try(confint(absDT, fitDT, nboot = 10),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})

absDF <- absoluteRisk(fitDF, time = 1)
absDT <- absoluteRisk(fitDF, time = 1)

test_that("no error in confint without newdata", {
    foo1 <- try(confint(absDF, fitDF, nboot = 10),
                silent = TRUE)
    foo2 <- try(confint(absDT, fitDT, nboot = 10),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})

absDF <- absoluteRisk(fitDF, time = 1, newdata = "typical")
absDT <- absoluteRisk(fitDF, time = 1, newdata = "typical")

test_that("no error in confint with typical profile", {
    foo1 <- try(confint(absDF, fitDF, nboot = 10),
                silent = TRUE)
    foo2 <- try(confint(absDT, fitDT, nboot = 10),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})

absDF <- absoluteRisk(fitDF, time = c(1, 2), newdata = DF[1, ])
absDT <- absoluteRisk(fitDF, time = c(1, 2), newdata = DT[1, ])

test_that("no error in confint with one covariate profile + 2 time points", {
    foo1 <- try(confint(absDF, fitDF, nboot = 10),
                silent = TRUE)
    foo2 <- try(confint(absDT, fitDT, nboot = 10),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})

absDF <- absoluteRisk(fitDF, time = c(1, 2), newdata = DF[c(1, n+1), ])
absDT <- absoluteRisk(fitDF, time = c(1, 2), newdata = DT[c(1, n+1), ])

test_that("no error in confint with two covariate profiles + 2 time points", {
    foo1 <- try(confint(absDF, fitDF, nboot = 10),
                silent = TRUE)
    foo2 <- try(confint(absDT, fitDT, nboot = 10),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})

# With splines
testthat::skip_if_not_installed("splines")

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
                         rep(1, n)),
                 "X2" = rnorm(2*n))

library(splines)
fit_obj <- fitSmoothHazard(event ~ Z*bs(X2), data = DF,
                           time = "ftime", ratio = 10)
times_for_absrisk <- seq(min(DF$ftime), max(DF$ftime), length.out = 50)

abs_obj <- absoluteRisk(fit_obj, time = times_for_absrisk,
                        newdata = DF[1, ],
                        method = "numerical")

test_that("no error in confint when using splines", {
    foo1 <- try(confint(abs_obj, fit_obj, nboot = 10),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
})

# Extrapolate
times_for_absrisk <- seq(min(DF$ftime), 2*max(DF$ftime), length.out = 50)
abs_obj <- absoluteRisk(fit_obj, time = times_for_absrisk,
                        newdata = DF[1, ],
                        method = "numerical")

test_that("no error in confint when extrapolating", {
    foo1 <- try(confint(abs_obj, fit_obj, nboot = 10),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
})
