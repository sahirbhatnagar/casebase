context("Absolute risk-Comp risk")

n = 100; alpha = 0.05

lambda10 <- 1
lambda20 <- 2
lambda11 <- 4
lambda21 <- 5

lambda_t0 <- lambda10 + lambda20
lambda_t1 <- lambda11 + lambda21

times <- c(rexp(n = n, rate = lambda_t0),
           rexp(n = n, rate = lambda_t1))
event <- c(rbinom(n, 1, prob = lambda10/lambda_t0),
           rbinom(n, 1, prob = lambda11/lambda_t1)) + 1
censor <- rexp(n = 2*n, rate = -log(alpha))

times_c <- pmin(times, censor)
event_c <- event * (times < censor)

DT <- data.table("ftime" = times_c,
                 "event" = event_c,
                 "Z" = c(rep(0,n), rep(1,n)))

DF <- data.frame("ftime" = times_c,
                 "event" = event_c,
                 "Z" = c(rep(0,n), rep(1,n)))

fitDF <- fitSmoothHazard(event ~ Z, data = DF, time = "ftime")
fitDT <- fitSmoothHazard(event ~ Z, data = DT, time = "ftime")

test_that("no error in absolute risk with data frames", {
    foo1 <- try(absoluteRisk(fitDF, time = 1, newdata = DF[1,],
                             method = "montecarlo"),
                silent = TRUE)
    foo2 <- try(absoluteRisk(fitDF, time = 1, newdata = DF[1,],
                             method = "numerical"),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})

test_that("no error in absolute risk with data tables", {
    foo1 <- try(absoluteRisk(fitDT, time = 1, newdata = DT[1,],
                             method = "montecarlo"),
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
    foo1 <- try(absoluteRisk(fitDF, time = 1, newdata = newDF,
                             method = "montecarlo"),
                silent = TRUE)
    foo2 <- try(absoluteRisk(fitDF, time = 1, newdata = newDF,
                             method = "numerical"),
                silent = TRUE)

    expect_false(inherits(foo1, "try-error"))
    expect_false(inherits(foo2, "try-error"))
})

test_that("no error in absolute risk with data tables - new data", {
    foo1 <- try(absoluteRisk(fitDT, time = 1, newdata = newDT,
                             method = "montecarlo"),
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

# Competing Risks and glmnet
extra_vars <- matrix(rnorm(10 * n), ncol = 10)
DF_ext <- cbind(DF, as.data.frame(extra_vars))
DT_ext <- cbind(DT, as.data.table(extra_vars))
formula_glmnet <- formula(paste(c("event ~ ftime", "Z",
                                  paste0("V", 1:10)),
                                collapse = " + "))

test_that("no error in fitting glmnet", {
    fitDF_glmnet <- try(fitSmoothHazard(formula_glmnet, data = DF_ext, time = "ftime",
                                        family = "glmnet", ratio = 10,
                                        lambda = c(0, 0.5)),
                        silent = TRUE)
    fitDT_glmnet <- try(fitSmoothHazard(formula_glmnet, data = DT_ext, time = "ftime",
                                        family = "glmnet", ratio = 10,
                                        lambda = c(0, 0.5)),
                        silent = TRUE)
    expect_false(inherits(fitDF_glmnet, "try-error"))
    expect_false(inherits(fitDT_glmnet, "try-error"))
})

y <- cbind(time = DF_ext$ftime,
           status = DF_ext$event)
x <- subset(DF_ext, select = !(colnames(DF_ext) %in% c("ftime", "event")))
x <- as.matrix(x)

test_that("no error in fitting fitSmoothHazard.fit", {
    fit_glmnet <- try(fitSmoothHazard.fit(x, y, time = "time", event = "status",
                                          family = "glmnet", ratio = 10,
                                          lambda = c(0, 0.5)),
                      silent = TRUE)

    expect_false(inherits(fit_glmnet, "try-error"))
})
