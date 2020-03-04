context("Fitting")

n = 100; alpha = 0.05

lambda_t0 <- 1
lambda_t1 <- 3

times <- c(rexp(n = n, rate = lambda_t0),
           rexp(n = n, rate = lambda_t1))
censor <- rexp(n = 2*n, rate = -log(alpha))

times_c <- pmin(times, censor)
event_c <- 1 * (times < censor)

DF <- data.frame("ftime" = times_c,
                 "event" = event_c,
                 "Z" = c(rep(0,n), rep(1,n)))
DT <- data.table("ftime" = times_c,
                 "event" = event_c,
                 "Z" = c(rep(0,n), rep(1,n)))

test_that("no error in fitting with data.frame and data.table", {
    fitDF <- try(fitSmoothHazard(event ~ Z, data = DF, time = "ftime"),
                silent = TRUE)
    fitDT <- try(fitSmoothHazard(event ~ Z, data = DT, time = "ftime"),
                 silent = TRUE)

    expect_false(inherits(fitDF, "try-error"))
    expect_false(inherits(fitDT, "try-error"))
})

# non-glm methods
extra_vars <- matrix(rnorm(10 * n), ncol = 10)
DF_ext <- cbind(DF, as.data.frame(extra_vars))
DT_ext <- cbind(DT, as.data.table(extra_vars))
formula_glmnet <- formula(paste(c("event ~ ftime", "Z",
                                  paste0("V", 1:10)),
                                collapse = " + "))

test_that("no error in fitting glmnet", {
    fitDF <- try(fitSmoothHazard(formula_glmnet, data = DF_ext, time = "ftime", family = "glmnet"),
                 silent = TRUE)
    fitDT <- try(fitSmoothHazard(formula_glmnet, data = DT_ext, time = "ftime", family = "glmnet"),
                 silent = TRUE)

    expect_false(inherits(fitDF, "try-error"))
    expect_false(inherits(fitDT, "try-error"))
})

formula_gam <- formula(paste(c("event ~ s(ftime)", "Z",
                                paste0("V", 1:10)),
                              collapse = " + "))

test_that("no error in fitting gam", {
    fitDF <- try(fitSmoothHazard(formula_gam, data = DF_ext, time = "ftime", family = "gam"),
                 silent = TRUE)
    fitDT <- try(fitSmoothHazard(formula_gam, data = DT_ext, time = "ftime", family = "gam"),
                 silent = TRUE)

    expect_false(inherits(fitDF, "try-error"))
    expect_false(inherits(fitDT, "try-error"))
})

formula_gbm <- formula(paste(c("event ~ ftime", "Z",
                               paste0("V", 1:10)),
                             collapse = " + "))
test_that("no error in fitting gbm", {
    fitDF <- try(fitSmoothHazard(formula_gbm, data = DF_ext, time = "ftime", family = "gbm"),
                 silent = TRUE)
    fitDT <- try(fitSmoothHazard(formula_gbm, data = DT_ext, time = "ftime", family = "gbm"),
                 silent = TRUE)

    expect_false(inherits(fitDF, "try-error"))
    expect_false(inherits(fitDT, "try-error"))
})

test_that("allow dot notation in formula", {
    try(model <- fitSmoothHazard(DeadOfPrCa ~ ., data = ERSPC, time='Follow.Up.Time', ratio = 100),
        silent = TRUE)

    expect_false(inherits(model, "try-error"))
})

test_that("sampling first and then fitting", {
    data_cb <- sampleCaseBase(ERSPC, time='Follow.Up.Time', ratio = 10, event = "DeadOfPrCa")
    try(model <- fitSmoothHazard(DeadOfPrCa ~ ., data = data_cb),
        silent = TRUE)

    expect_false(inherits(model, "try-error"))
})
