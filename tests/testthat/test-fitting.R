context("Fitting")

n <- 100
alp <- 0.05
lambda_t0 <- 1
lambda_t1 <- 3

times <- c(rexp(n = n, rate = lambda_t0),
           rexp(n = n, rate = lambda_t1))
censor <- rexp(n = 2 * n, rate = -log(alp))

times_c <- pmin(times, censor)
event_c <- 1 * (times < censor)

DF <- data.frame("ftime" = times_c,
                 "event" = event_c,
                 "Z" = c(rep(0, n),
                         rep(1, n)))
DT <- data.table("ftime" = times_c,
                 "event" = event_c,
                 "Z" = c(rep(0, n),
                         rep(1, n)))

test_that("no error in fitting with data.frame and data.table", {
    fitDF <- try(fitSmoothHazard(event ~ Z, data = DF, time = "ftime"),
                silent = TRUE)
    fitDT <- try(fitSmoothHazard(event ~ Z, data = DT, time = "ftime"),
                 silent = TRUE)

    expect_false(inherits(fitDF, "try-error"))
    expect_false(inherits(fitDT, "try-error"))
})

test_that("allow dot notation in formula", {
    try(model <- fitSmoothHazard(DeadOfPrCa ~ ., data = ERSPC,
                                 time = "Follow.Up.Time", ratio = 100),
        silent = TRUE)

    expect_false(inherits(model, "try-error"))
})

test_that("sampling first and then fitting", {
    data_cb <- sampleCaseBase(ERSPC, time = "Follow.Up.Time",
                              ratio = 10, event = "DeadOfPrCa")
    try(model <- fitSmoothHazard(DeadOfPrCa ~ ., data = data_cb),
        silent = TRUE)

    expect_false(inherits(model, "try-error"))
})

#####################
# Formula parsing----
form <- formula(event ~ exposure + time)
form_bs <- formula(event ~ exposure + bs(time))
form_log <- formula(event ~ exposure + log(time))
form_int <- formula(event ~ exposure * time)
form_nested <- formula(cens ~ horTh * nsx(log(time), df = 3) + age * time)

form_bs_extra <- formula(event ~ exposure + bs(time, df = 3))
form_bs_named <- formula(event ~ exposure + bs(x = time))

test_that("detecting non-linear functions of time", {
    expect_false(detect_nonlinear_time(form, "time"))
    expect_true(detect_nonlinear_time(form_bs, "time"))
    expect_true(detect_nonlinear_time(form_log, "time"))
    expect_false(detect_nonlinear_time(form_int, "time"))
    expect_true(detect_nonlinear_time(form_nested, "time"))

    expect_true(detect_nonlinear_time(form_bs_extra, "time"))
    expect_true(detect_nonlinear_time(form_bs_named, "time"))
})

wrong <- formula(event ~ exposure + time + wrongtime)
wrong_bs <- formula(event ~ exposure + time + bs(wrongtime))
wrong_log <- formula(event ~ exposure + time + log(wrongtime))
wrong_int <- formula(event ~ exposure * wrongtime + time)

wrong2 <- formula(event ~ exposure + time + timewrong)
wrong2_bs <- formula(event ~ exposure + time + bs(timewrong))
wrong2_log <- formula(event ~ exposure + time + log(timewrong))
wrong2_int <- formula(event ~ exposure * timewrong + time)

test_that("Making sure we don't pick up anything that looks like time", {
    expect_false(detect_nonlinear_time(wrong, "time"))
    expect_false(detect_nonlinear_time(wrong_bs, "time"))
    expect_false(detect_nonlinear_time(wrong_log, "time"))
    expect_false(detect_nonlinear_time(wrong_int, "time"))
    expect_false(detect_nonlinear_time(wrong2, "time"))
    expect_false(detect_nonlinear_time(wrong2_bs, "time"))
    expect_false(detect_nonlinear_time(wrong2_log, "time"))
    expect_false(detect_nonlinear_time(wrong2_int, "time"))
})

test_that("detecting interactions with time", {
    expect_false(detect_interaction(form))
    expect_false(detect_interaction(form_bs))
    expect_false(detect_interaction(form_log))
    expect_true(detect_interaction(form_int))
})
