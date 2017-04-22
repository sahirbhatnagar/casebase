context("Absolute risk")

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
