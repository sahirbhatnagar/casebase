context("Absolute risk")

# Handling warning messages coming from montecarlo integration
handler_validmc <- function(msg) {
    if (any(grepl("out of range", msg))) invokeRestart("muffleWarning")
}

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

fitDF <- fitSmoothHazard(event ~ Z, data = DF, time = "ftime", ratio = 10)
fitDT <- fitSmoothHazard(event ~ Z, data = DT, time = "ftime", ratio = 10)

test_that("no error in absolute risk with data frames", {
    foo1 <- try(withCallingHandlers(absoluteRisk(fitDF, time = 1, newdata = DF[1,],
                                                 method = "montecarlo"),
                                    warning = handler_validmc),
                silent = TRUE)
    foo2 <- try(absoluteRisk(fitDF, time = 1, newdata = DF[1,],
                             method = "numerical"),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})

test_that("no error in absolute risk with data tables", {
    foo1 <- try(withCallingHandlers(absoluteRisk(fitDT, time = 1, newdata = DT[1,],
                                                 method = "montecarlo"),
                                    warning = handler_validmc),
                silent = TRUE)
    foo2 <- try(absoluteRisk(fitDT, time = 1, newdata = DT[1,],
                             method = "numerical"),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})

# Using new data
newDT <- data.table("Z" = c(0,1))
newDF <- data.frame("Z" = c(0,1))

test_that("no error in absolute risk with data frames - new data", {
    foo1 <- try(withCallingHandlers(absoluteRisk(fitDF, time = 1, newdata = newDF,
                             method = "montecarlo"),
                             warning = handler_validmc),
                silent = TRUE)
    foo2 <- try(absoluteRisk(fitDF, time = 1, newdata = newDF,
                             method = "numerical"),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})

test_that("no error in absolute risk with data tables - new data", {
    foo1 <- try(withCallingHandlers(absoluteRisk(fitDT, time = 1, newdata = newDT,
                             method = "montecarlo"),
                             warning = handler_validmc),
                silent = TRUE)
    foo2 <- try(absoluteRisk(fitDT, time = 1, newdata = newDT,
                             method = "numerical"),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})


# Make sure we get probabilities
test_that("should output probabilities with data frames", {
    absRiskMC <- absoluteRisk(fitDF, time = 1, newdata = DF[1,],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDF, time = 1, newdata = DF[1,],
                              method = "numerical")

    expect_true(all(absRiskMC >= 0))
    expect_true(all(absRiskNI >= 0))
    expect_true(all(absRiskMC <= 1))
    expect_true(all(absRiskNI <= 1))
})

test_that("should output probabilities with data tables", {
    absRiskMC <- absoluteRisk(fitDT, time = 1, newdata = DT[1,],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDT, time = 1, newdata = DT[1,],
                              method = "numerical")

    expect_true(all(absRiskMC >= 0))
    expect_true(all(absRiskNI >= 0))
    expect_true(all(absRiskMC <= 1))
    expect_true(all(absRiskNI <= 1))
})

test_that("should output probabilities with data frames - two time points", {
    absRiskMC <- absoluteRisk(fitDF, time = c(0.5, 1), newdata = DF[1,],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDF, time = c(0.5, 1), newdata = DF[1,],
                              method = "numerical")

    expect_true(all(absRiskMC[,-1] >= 0))
    expect_true(all(absRiskNI[,-1] >= 0))
    expect_true(all(absRiskMC[,-1] <= 1))
    expect_true(all(absRiskNI[,-1] <= 1))
})

test_that("should output probabilities with data tables - two time points", {
    absRiskMC <- absoluteRisk(fitDT, time = c(0.5, 1), newdata = DT[1,],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDT, time = c(0.5, 1), newdata = DT[1,],
                              method = "numerical")

    expect_true(all(absRiskMC[,-1] >= 0))
    expect_true(all(absRiskNI[,-1] >= 0))
    expect_true(all(absRiskMC[,-1] <= 1))
    expect_true(all(absRiskNI[,-1] <= 1))
})

test_that("should output probabilities with data frames - two covariate profile", {
    absRiskMC <- absoluteRisk(fitDF, time = 1, newdata = DF[c(1, n + 1),],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDF, time = 1, newdata = DF[c(1, n + 1),],
                              method = "numerical")

    expect_true(all(absRiskMC >= 0))
    expect_true(all(absRiskNI >= 0))
    expect_true(all(absRiskMC <= 1))
    expect_true(all(absRiskNI <= 1))
})

test_that("should output probabilities with data tables - two covariate profile", {
    absRiskMC <- absoluteRisk(fitDT, time = 1, newdata = DT[c(1, n + 1),],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDT, time = 1, newdata = DT[c(1, n + 1),],
                              method = "numerical")

    expect_true(all(absRiskMC >= 0))
    expect_true(all(absRiskNI >= 0))
    expect_true(all(absRiskMC <= 1))
    expect_true(all(absRiskNI <= 1))
})

# Absolute risk at time = 0
test_that("should give probability 0 at time 0 with data frames", {
    absRiskMC <- absoluteRisk(fitDF, time = 0, newdata = newDF,
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDF, time = 0, newdata = newDF,
                              method = "numerical")

    expect_true(all.equal(absRiskMC,
                          matrix(0, ncol = 2, nrow = 1),
                          check.attributes = FALSE))
    expect_true(all.equal(absRiskNI,
                          matrix(0, ncol = 2, nrow = 1),
                          check.attributes = FALSE))
})

test_that("should give probability 0 at time 0 with data tables", {
    absRiskMC <- absoluteRisk(fitDT, time = 0, newdata = newDT,
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDT, time = 0, newdata = newDT,
                              method = "numerical")

    expect_true(all.equal(absRiskMC,
                          matrix(0, ncol = 2, nrow = 1),
                          check.attributes = FALSE))
    expect_true(all.equal(absRiskNI,
                          matrix(0, ncol = 2, nrow = 1),
                          check.attributes = FALSE))
})

test_that("should compute risk when time and newdata aren't provided", {
    absRiskDF <- absoluteRisk(fitDF)
    absRiskDT <- absoluteRisk(fitDT)

    expect_true("risk" %in% names(absRiskDF))
    expect_true("risk" %in% names(absRiskDT))
})

# non-glm methods----
fitDF_gam <- fitSmoothHazard(event ~ s(ftime) + Z, data = DF, time = "ftime", family = "gam", ratio = 10)
fitDT_gam <- fitSmoothHazard(event ~ s(ftime) + Z, data = DT, time = "ftime", family = "gam", ratio = 10)

test_that("no error in fitting gam", {
    riskDF <- try(absoluteRisk(fitDF_gam, time = 0.5, newdata = newDF),
                  silent = TRUE)
    riskDT <- try(absoluteRisk(fitDT_gam, time = 0.5, newdata = newDT),
                  silent = TRUE)

    expect_false(inherits(riskDF, "try-error"))
    expect_false(inherits(riskDT, "try-error"))
})

test_that("should compute risk when time and newdata aren't provided", {
    absRiskDF_gam <- absoluteRisk(fitDF_gam)
    absRiskDT_gam <- absoluteRisk(fitDT_gam)

    expect_true("risk" %in% names(absRiskDF_gam))
    expect_true("risk" %in% names(absRiskDT_gam))
})

fitDF_gbm <- fitSmoothHazard(event ~ ftime + Z, data = DF, time = "ftime", family = "gbm", ratio = 10)
fitDT_gbm <- fitSmoothHazard(event ~ ftime + Z, data = DT, time = "ftime", family = "gbm", ratio = 10)

test_that("no error in fitting gbm", {
    riskDF <- try(absoluteRisk(fitDF_gbm, time = 0.5, newdata = newDF, n.trees = 100, nsamp = 500),
                 silent = TRUE)
    riskDT <- try(absoluteRisk(fitDT_gbm, time = 0.5, newdata = newDT, n.trees = 100, nsamp = 500),
                 silent = TRUE)
    riskDF_mc <- try(absoluteRisk(fitDF_gbm, time = 0.5, newdata = newDF, n.trees = 100, nsamp = 10,
                                  method = "montecarlo"),
                     silent = TRUE)
    riskDT_mc <- try(absoluteRisk(fitDT_gbm, time = 0.5, newdata = newDT, n.trees = 100, nsamp = 10,
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
    fitDF_gbm_red$originalData <- fitDF_gbm$originalData[c(1:5, 101:105),]
    absRiskDF_gbm <- absoluteRisk(fitDF_gbm_red, n.trees = 100, nsamp = 500)

    fitDT_gbm_red <- fitDT_gbm
    fitDT_gbm_red$originalData <- fitDT_gbm$originalData[c(1:5, 101:105),]
    absRiskDT_gbm <- absoluteRisk(fitDT_gbm_red, n.trees = 100, nsamp = 500)

    expect_true("risk" %in% names(absRiskDF_gbm))
    expect_true("risk" %in% names(absRiskDT_gbm))
})

extra_vars <- matrix(rnorm(10 * n), ncol = 10)
DF_ext <- cbind(DF, as.data.frame(extra_vars))
DT_ext <- cbind(DT, as.data.table(extra_vars))
formula_glmnet <- formula(paste(c("event ~ ftime", "Z",
                                  paste0("V", 1:10)),
                                collapse = " + "))
fitDF_glmnet <- fitSmoothHazard(formula_glmnet, data = DF_ext, time = "ftime", family = "glmnet", ratio = 10)
fitDT_glmnet <- fitSmoothHazard(formula_glmnet, data = DT_ext, time = "ftime", family = "glmnet", ratio = 10)

extra_vars_new <- matrix(rnorm(10 * 2), ncol = 10)
colnames(extra_vars_new) <- paste0("V", 1:10)
newDF_ext <- cbind(newDF, extra_vars_new)
newDT_ext <- cbind(newDT, extra_vars_new)

test_that("no error in fitting glmnet", {
    riskDF <- try(absoluteRisk(fitDF_glmnet, time = 0.5, newdata = newDF_ext),
                  silent = TRUE)
    riskDT <- try(absoluteRisk(fitDT_glmnet, time = 0.5, newdata = newDT_ext),
                  silent = TRUE)
    riskDF_mc <- try(absoluteRisk(fitDF_glmnet, time = 0.5, newdata = newDF_ext, nsamp = 10,
                                  method = "montecarlo"),
                     silent = TRUE)
    riskDT_mc <- try(absoluteRisk(fitDT_glmnet, time = 0.5, newdata = newDT_ext, nsamp = 10,
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
    riskDF_gam <- absoluteRisk(fitDF_gam, time = 0.5, newdata = newDF, family = "gam")
    riskDT_gam <- absoluteRisk(fitDT_gam, time = 0.5, newdata = newDT, family = "gam")

    riskDF_gbm <- absoluteRisk(fitDF_gbm, time = 0.5, newdata = newDF,
                               family = "gbm", n.trees = 100, nsamp = 500)
    riskDT_gbm <- absoluteRisk(fitDT_gbm, time = 0.5, newdata = newDT,
                               family = "gbm", n.trees = 100, nsamp = 500)

    riskDF_glmnet <- absoluteRisk(fitDF_glmnet, time = 0.5, newdata = newDF_ext, family = "glmnet")
    riskDT_glmnet <- absoluteRisk(fitDT_glmnet, time = 0.5, newdata = newDT_ext, family = "glmnet")

    expect_true(all(riskDF_gam >= 0))
    expect_true(all(riskDT_gam >= 0))
    expect_true(all(riskDF_gam <= 1))
    expect_true(all(riskDT_gam <= 1))

    expect_true(all(riskDF_gbm >= 0))
    expect_true(all(riskDT_gbm >= 0))
    expect_true(all(riskDF_gbm <= 1))
    expect_true(all(riskDT_gbm <= 1))

    expect_true(all(riskDF_glmnet >= 0))
    expect_true(all(riskDT_glmnet >= 0))
    expect_true(all(riskDF_glmnet <= 1))
    expect_true(all(riskDT_glmnet <= 1))
})

# test_that("should compute risk when time and newdata aren't provided", {
#     absRiskDF <- absoluteRisk(fitDF)
#     absRiskDT <- absoluteRisk(fitDT)
#
#     expect_true("risk" %in% names(absRiskDF))
#     expect_true("risk" %in% names(absRiskDT))
# })
