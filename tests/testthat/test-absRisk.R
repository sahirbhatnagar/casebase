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
    foo <- try(absoluteRisk(fitDF, time = 1, newdata = DF[1,]))

    expect_false(inherits(foo, "try-error"))
})

test_that("no error in absolute risk with data tables", {
    foo <- try(absoluteRisk(fitDT, time = 1, newdata = DT[1,]))

    expect_false(inherits(foo, "try-error"))
})

# Using new data
newDT <- data.table("Z" = c(0,1))
newDF <- data.frame("Z" = c(0,1))

test_that("no error in absolute risk with data frames - new data", {
    foo <- try(absoluteRisk(fitDF, time = 1, newdata = newDF))

    expect_false(inherits(foo, "try-error"))
})

test_that("no error in absolute risk with data tables - new data", {
    foo <- try(absoluteRisk(fitDT, time = 1, newdata = newDT))

    expect_false(inherits(foo, "try-error"))
})

# Make sure we get probabilities
test_that("should output probabilities with data frames", {
    absRiskMC <- absoluteRisk(fitDF, time = 1, newdata = DF[1,],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDF, time = 1, newdata = DF[1,],
                              method = "quadrature")

    expect_true(all(absRiskMC >= 0))
    expect_true(all(absRiskNI >= 0))
    expect_true(all(absRiskMC <= 1))
    expect_true(all(absRiskNI <= 1))
})

test_that("should output probabilities with data tables", {
    absRiskMC <- absoluteRisk(fitDT, time = 1, newdata = DT[1,],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDT, time = 1, newdata = DT[1,],
                              method = "quadrature")

    expect_true(all(absRiskMC >= 0))
    expect_true(all(absRiskNI >= 0))
    expect_true(all(absRiskMC <= 1))
    expect_true(all(absRiskNI <= 1))
})

test_that("should output probabilities with data frames - two time points", {
    absRiskMC <- absoluteRisk(fitDF, time = c(0.5, 1), newdata = DF[1,],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDF, time = c(0.5, 1), newdata = DF[1,],
                              method = "quadrature")

    expect_true(all(absRiskMC >= 0))
    expect_true(all(absRiskNI >= 0))
    expect_true(all(absRiskMC <= 1))
    expect_true(all(absRiskNI <= 1))
})

test_that("should output probabilities with data tables - two time points", {
    absRiskMC <- absoluteRisk(fitDT, time = c(0.5, 1), newdata = DT[1,],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDT, time = c(0.5, 1), newdata = DT[1,],
                              method = "quadrature")

    expect_true(all(absRiskMC >= 0))
    expect_true(all(absRiskNI >= 0))
    expect_true(all(absRiskMC <= 1))
    expect_true(all(absRiskNI <= 1))
})

test_that("should output probabilities with data frames - two covariate profile", {
    absRiskMC <- absoluteRisk(fitDF, time = 1, newdata = DF[c(1, n + 1),],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDF, time = 1, newdata = DF[c(1, n + 1),],
                              method = "quadrature")

    expect_true(all(absRiskMC >= 0))
    expect_true(all(absRiskNI >= 0))
    expect_true(all(absRiskMC <= 1))
    expect_true(all(absRiskNI <= 1))
})

test_that("should output probabilities with data tables - two covariate profile", {
    absRiskMC <- absoluteRisk(fitDT, time = 1, newdata = DT[c(1, n + 1),],
                              method = "montecarlo")
    absRiskNI <- absoluteRisk(fitDT, time = 1, newdata = DT[c(1, n + 1),],
                              method = "quadrature")

    expect_true(all(absRiskMC >= 0))
    expect_true(all(absRiskNI >= 0))
    expect_true(all(absRiskMC <= 1))
    expect_true(all(absRiskNI <= 1))
})