library(casebase)
context("Absolute risk")

test_that("no error in absolute risk", {
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
    DT <- data.frame("ftime" = times_c,
                     "event" = event_c,
                     "Z" = c(rep(0,n), rep(1,n)))

    fit <- fitSmoothHazard(event ~ Z, data=DT, time="ftime")
    foo <- try(absoluteRisk(fit, time = 1, newdata = DT[1,]))

    expect_false(inherits(foo, "try-error"))
})