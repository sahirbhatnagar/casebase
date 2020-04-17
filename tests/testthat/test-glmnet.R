context("glmnet")

# Skip tests if gbm is not installed
testthat::skip_if_not_installed("glmnet")

# Create data----
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

extra_vars <- matrix(rnorm(10 * n), ncol = 10)
DF_ext <- cbind(DF, as.data.frame(extra_vars))
DT_ext <- cbind(DT, as.data.table(extra_vars))

formula_glmnet <- formula(paste(c("event ~ ftime", "Z",
                                  paste0("V", 1:10)),
                                collapse = " + "))

# Fitting----
test_that("no error in fitting glmnet", {
    fitDF <- try(fitSmoothHazard(formula_glmnet, data = DF_ext, time = "ftime", family = "glmnet"),
                 silent = TRUE)
    fitDT <- try(fitSmoothHazard(formula_glmnet, data = DT_ext, time = "ftime", family = "glmnet"),
                 silent = TRUE)

    expect_false(inherits(fitDF, "try-error"))
    expect_false(inherits(fitDT, "try-error"))
})

# Absolute risk----
fitDF_glmnet <- fitSmoothHazard(formula_glmnet, data = DF_ext, time = "ftime", family = "glmnet", ratio = 10)
fitDT_glmnet <- fitSmoothHazard(formula_glmnet, data = DT_ext, time = "ftime", family = "glmnet", ratio = 10)

newDT <- data.table("Z" = c(0,1))
newDF <- data.frame("Z" = c(0,1))

extra_vars_new <- matrix(rnorm(10 * 2), ncol = 10)
colnames(extra_vars_new) <- paste0("V", 1:10)
newDF_ext <- cbind(newDF, extra_vars_new)
newDT_ext <- cbind(newDT, extra_vars_new)


test_that("no error in fitting glmnet", {
    riskDF <- try(absoluteRisk(fitDF_glmnet, time = 0.5, newdata = newDF_ext),
                  silent = TRUE)
    riskDT <- try(absoluteRisk(fitDT_glmnet, time = 0.5, newdata = newDT_ext),
                  silent = TRUE)
    riskDF_mc <- try(absoluteRisk(fitDF_glmnet, time = 0.5, newdata = newDF_ext,
                                  nsamp = 10,
                                  method = "montecarlo"),
                     silent = TRUE)
    riskDT_mc <- try(absoluteRisk(fitDT_glmnet, time = 0.5, newdata = newDT_ext,
                                  nsamp = 10,
                                  method = "montecarlo"),
                     silent = TRUE)

    expect_false(inherits(riskDF, "try-error"))
    expect_false(inherits(riskDT, "try-error"))
    expect_false(inherits(riskDF_mc, "try-error"))
    expect_false(inherits(riskDT_mc, "try-error"))
})

test_that("no error in using custom lambda in glmnet", {
    riskDF <- try(absoluteRisk(fitDF_glmnet, time = 0.5, newdata = newDF_ext, s = 0.1),
                  silent = TRUE)
    riskDT <- try(absoluteRisk(fitDT_glmnet, time = 0.5, newdata = newDT_ext, s = 0.1),
                  silent = TRUE)

    expect_false(inherits(riskDF, "try-error"))
    expect_false(inherits(riskDT, "try-error"))
})

test_that("output probabilities", {
    riskDF_glmnet <- absoluteRisk(fitDF_glmnet, time = 0.5, newdata = newDF_ext, family = "glmnet")
    riskDT_glmnet <- absoluteRisk(fitDT_glmnet, time = 0.5, newdata = newDT_ext, family = "glmnet")

    expect_true(all(riskDF_glmnet >= 0))
    expect_true(all(riskDT_glmnet >= 0))
    expect_true(all(riskDF_glmnet <= 1))
    expect_true(all(riskDT_glmnet <= 1))
})

# Matrix interface----
N <- 1000; p <- 30
nzc <- p/3
x <- matrix(rnorm(N*p),N,p)
dimnames(x)[[2]] <- paste0("x",1:p)
beta <- rnorm(nzc)
fx <- x[,seq(nzc)] %*% beta/3
hx <- exp(fx)
ty <- rexp(N,hx)
tcens <- rbinom(n = N,
                prob = 0.3,
                size = 1) # censoring indicator
y <- cbind(time = ty, status = 1 - tcens) # y=Surv(ty,1-tcens) with library(survival)

test_that("no error in fitting fitSmoothHazard.fit", {
    fit_glmnet <- try(fitSmoothHazard.fit(x, y, time = "time", event = "status",
                                          family = "glmnet", ratio = 10,
                                          lambda = c(0, 0.5)),
                      silent = TRUE)

    expect_false(inherits(fit_glmnet, "try-error"))
})

test_that("no error in using nonlinear functions of time", {
    skip_if_not_installed("splines")
    library(splines)
    fit_glmnet <- try(fitSmoothHazard.fit(x, y, formula_time = ~ log(time),
                                          time = "time", event = "status",
                                          family = "glmnet", ratio = 10,
                                          lambda = c(0, 0.5)),
                      silent = TRUE)

    fit_glmnet_splines <- try(fitSmoothHazard.fit(x, y, formula_time = ~ bs(time),
                                                  time = "time", event = "status",
                                                  family = "glmnet", ratio = 10,
                                                  lambda = c(0, 0.5)),
                              silent = TRUE)

    expect_false(inherits(fit_glmnet, "try-error"))
    expect_false(inherits(fit_glmnet_splines, "try-error"))
})

# Absolute risk with matrix interface----
# Only family="glmnet" has been implemented
fit_glmnet <- fitSmoothHazard.fit(x, y, time = "time", event = "status",
                                  family = "glmnet", ratio = 10,
                                  lambda = c(0, 0.5))
fit_glmnet_log <- fitSmoothHazard.fit(x, y, formula_time = ~ log(time),
                                      time = "time", event = "status",
                                      family = "glmnet", ratio = 10,
                                      lambda = c(0, 0.5))
new_x <- x[1:10, ]
risk <- try(absoluteRisk(fit_glmnet, time = 1,
                         newdata = new_x, nsamp = 100),
            silent = TRUE)
risk_log <- try(absoluteRisk(fit_glmnet_log, time = 1,
                             newdata = new_x, nsamp = 100),
                silent = TRUE)


test_that("no error in absoluteRisk with glmnet", {
    expect_false(inherits(risk, "try-error"))
    expect_false(inherits(risk_log, "try-error"))
})

test_that("we get probabilities", {
    expect_true(all(risk >= 0))
    expect_true(all(risk <= 1))
    expect_true(all(risk_log >= 0))
    expect_true(all(risk_log <= 1))
})

# Test absoluteRisk--two time points
risk <- try(absoluteRisk(fit_glmnet, time = c(1,2),
                         newdata = new_x, nsamp = 100),
            silent = TRUE)
risk_log <- try(absoluteRisk(fit_glmnet_log, time = c(1,2),
                             newdata = new_x, nsamp = 100),
                silent = TRUE)

test_that("no error in absoluteRisk with glmnet", {

    expect_false(inherits(risk, "try-error"))
    expect_false(inherits(risk_log, "try-error"))
})

test_that("we get probabilities", {
    expect_true(all(risk[,-1] >= 0))
    expect_true(all(risk[,-1] <= 1))
    expect_true(all(risk_log[,-1] >= 0))
    expect_true(all(risk_log[,-1] <= 1))
})
