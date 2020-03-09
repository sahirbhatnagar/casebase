context("Matrix interface")
library(splines)

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
    fit_glm <- try(fitSmoothHazard.fit(x, y, time = "time", event = "status", ratio = 10),
                   silent = TRUE)
    fit_glmnet <- try(fitSmoothHazard.fit(x, y, time = "time", event = "status",
                                          family = "glmnet", ratio = 10,
                                          lambda = c(0, 0.5)),
                      silent = TRUE)
    fit_gbm <- try(fitSmoothHazard.fit(x, y, time = "time", event = "status",
                                       family = "gbm", ratio = 10),
                   silent = TRUE)

    expect_false(inherits(fit_glm, "try-error"))
    expect_false(inherits(fit_glmnet, "try-error"))
    expect_false(inherits(fit_gbm, "try-error"))
})

test_that("no error in using nonlinear functions of time", {
    fit_glm <- try(fitSmoothHazard.fit(x, y, formula_time = ~ log(time),
                                       time = "time", event = "status", ratio = 10),
                   silent = TRUE)
    fit_glmnet <- try(fitSmoothHazard.fit(x, y, formula_time = ~ log(time),
                                          time = "time", event = "status",
                                          family = "glmnet", ratio = 10,
                                          lambda = c(0, 0.5)),
                      silent = TRUE)
    # fit_gbm <- try(fitSmoothHazard.fit(x, y, formula_time = ~ log(time),
    #                                    time = "time", event = "status",
    #                                    family = "gbm", ratio = 10),
    #                silent = TRUE)

    fit_glm_splines <- try(fitSmoothHazard.fit(x, y, formula_time = ~ bs(time),
                                               time = "time", event = "status", ratio = 10),
                           silent = TRUE)
    fit_glmnet_splines <- try(fitSmoothHazard.fit(x, y, formula_time = ~ bs(time),
                                                  time = "time", event = "status",
                                                  family = "glmnet", ratio = 10,
                                                  lambda = c(0, 0.5)),
                              silent = TRUE)
    # fit_gbm_splines <- try(fitSmoothHazard.fit(x, y, formula_time = ~ bs(time),
    #                                            time = "time", event = "status",
    #                                            family = "gbm", ratio = 10),
    #                        silent = TRUE)

    expect_false(inherits(fit_glm, "try-error"))
    expect_false(inherits(fit_glmnet, "try-error"))
    # expect_false(inherits(fit_gbm, "try-error"))

    expect_false(inherits(fit_glm_splines, "try-error"))
    expect_false(inherits(fit_glmnet_splines, "try-error"))
    # expect_false(inherits(fit_gbm_splines, "try-error"))
})

test_that("warnings when using gbm witn non-linear functions of time or interactions", {
    expect_warning(fitSmoothHazard.fit(x, y, formula_time = ~ log(time),
                                       time = "time", event = "status",
                                       family = "gbm", ratio = 10),
                   regexp = "gbm may throw an error")
    expect_warning(fitSmoothHazard.fit(x, y, formula_time = ~ bs(time),
                                       time = "time", event = "status",
                                       family = "gbm", ratio = 10),
                   regexp = "gbm may throw an error")
})

fit_glmnet <- fitSmoothHazard.fit(x, y, time = "time", event = "status",
                                  family = "glmnet", ratio = 10,
                                  lambda = c(0, 0.5))
fit_glmnet_log <- fitSmoothHazard.fit(x, y, formula_time = ~ log(time),
                                      time = "time", event = "status",
                                      family = "glmnet", ratio = 10,
                                      lambda = c(0, 0.5))
# Test absoluteRisk----
# Only family="glmnet" has been implemented
new_x <- x[1:10, ]
risk <- try(absoluteRisk(fit_glmnet, time = 1,
                         newdata = new_x, nsamp = 100),
            silent = TRUE)
risk_log <- try(absoluteRisk(fit_glmnet_log, time = 1,
                             newdata = new_x, nsamp = 100),
                silent = TRUE)

test_that("error with glm.fit",{
    expect_error(absoluteRisk(fit_glm, time = 1,
                              newdata = new_x, nsamp = 100))
})

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
