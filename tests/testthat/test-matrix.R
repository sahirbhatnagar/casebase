context("Matrix interface")
set.seed(12345)

# CRAN skip atlas check fix
testthat::skip_if(grepl(pattern = "atlas", sessionInfo()$BLAS,
                        ignore.case = TRUE))

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
y <- cbind(time = ty, status = 1 - tcens)

test_that("no error in fitting fitSmoothHazard.fit", {
    fit_glm <- try(fitSmoothHazard.fit(x, y, time = "time", event = "status",
                                       ratio = 10),
                   silent = TRUE)

    expect_false(inherits(fit_glm, "try-error"))
})

test_that("no error in using nonlinear functions of time", {
    skip_if_not_installed("splines")
    library(splines)
    fit_glm <- try(fitSmoothHazard.fit(x, y, formula_time = ~ log(time),
                                       time = "time", event = "status",
                                       ratio = 10),
                   silent = TRUE)
    fit_glm_splines <- try(fitSmoothHazard.fit(x, y, formula_time = ~ bs(time),
                                               time = "time", event = "status",
                                               ratio = 10),
                           silent = TRUE)

    expect_false(inherits(fit_glm, "try-error"))
    expect_false(inherits(fit_glm_splines, "try-error"))
})

test_that("error with glm.fit", {
    expect_error(absoluteRisk(fit_glm, time = 1,
                              newdata = new_x, nsamp = 100))
})
