context("Fitting-Comp risk")
set.seed(12345)

# CRAN skip atlas check fix
testthat::skip_if(grepl(pattern = "atlas", sessionInfo()$BLAS,
                        ignore.case = TRUE))


n <- 100
alp <- 0.05
lambda10 <- 1
lambda20 <- 2
lambda11 <- 4
lambda21 <- 5

lambda_t0 <- lambda10 + lambda20
lambda_t1 <- lambda11 + lambda21

times <- c(rexp(n = n, rate = lambda_t0),
           rexp(n = n, rate = lambda_t1))
event <- c(rbinom(n, 1, prob = lambda10 / lambda_t0),
           rbinom(n, 1, prob = lambda11 / lambda_t1)) + 1
censor <- rexp(n = 2 * n, rate = -log(alp))

times_c <- pmin(times, censor)
event_c <- event * (times < censor)

DT <- data.table("ftime" = times_c,
                 "event" = event_c,
                 "Z" = c(rep(0, n),
                         rep(1, n)))

DF <- data.frame("ftime" = times_c,
                 "event" = event_c,
                 "Z" = c(rep(0, n),
                         rep(1, n)))

test_that("no error in fitting with data.frame and data.table", {
    fitDF <- try(fitSmoothHazard(event ~ Z, data = DF, time = "ftime"),
                 silent = TRUE)
    fitDT <- try(fitSmoothHazard(event ~ Z, data = DT, time = "ftime"),
                 silent = TRUE)

    expect_false(inherits(fitDF, "try-error"))
    expect_false(inherits(fitDF, "try-error"))
})

test_that("Errors in fitting when family is not glm or glmnet", {
    expect_error(fitSmoothHazard(event ~ Z, data = DF, time = "ftime",
                                 family = "gam"))
    expect_error(fitSmoothHazard(event ~ Z, data = DF, time = "ftime",
                                 family = "gbm"))
    expect_error(fitSmoothHazard(event ~ Z, data = DF, time = "ftime",
                                 family = "other_family"))
})

# Sample first then fit----
# See issue 149
dta <- data.frame(
    time = exp(rnorm(100)),
    event = sample(0:2, size = 100, replace = TRUE)
)

# create dataset with casebase samples
cbdata <- sampleCaseBase(data = dta, time = "time", event = "event",
                         comprisk = TRUE)

#try to fit the smooth hazard model for competing events
test_that("No error when sampling first then fitting", {
    fit <- try(fitSmoothHazard(event ~ time, data = cbdata, time = "time"),
               silent = TRUE)

    expect_false(inherits(fit, "try-error"))
})

