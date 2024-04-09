context("popTime methods")
set.seed(12345)

# CRAN skip atlas check fix
testthat::skip_if(grepl(pattern = "atlas", sessionInfo()$BLAS,
                        ignore.case = TRUE))

nobs <- 500

# simulation parameters
a1 <- 1.0
b1 <- 200
a2 <- 1.0
b2 <- 50
c1 <- 0.0
c2 <- 0.0

# end of study time
eost <- 10

# e event type 0-censored, 1-event of interest, 2-competing event
# t observed time/endpoint
# z is a binary covariate
DTsim <- data.table(ID = seq_len(nobs), z = rbinom(nobs, 1, 0.5))
data.table::setkey(DTsim, ID)
DTsim[, `:=`(event_time = rweibull(nobs, a1, b1 * exp(z * c1)),
             competing_time = rweibull(nobs, a2, b2 * exp(z * c2)),
             end_of_study_time = eost)]
DTsim[, `:=`(event = 1 * (event_time < competing_time) +
                2 * (event_time >= competing_time),
            time = pmin(event_time, competing_time))]
DTsim[time >= end_of_study_time, event := 0]
DTsim[time >= end_of_study_time, time := end_of_study_time]

DFsim <- data.frame("z" = DTsim[, z],
                    "event" = DTsim[, event],
                    "time" = DTsim[, time])

# data.frames input
out1 <- popTime(data = DFsim, time = "time", event = "event")
out2 <- popTime(data = DFsim, time = "time", event = "event", exposure = "z")

# data.table input
out3 <- popTime(data = DTsim, time = "time", event = "event")
out4 <- popTime(data = DTsim, time = "time", event = "event", exposure = "z")

test_that("expect data.table or data.frame in popTime with data.frame or data.table input", {

    expect_s3_class(out1, "data.frame")
    expect_s3_class(out1, "data.table")
    expect_s3_class(out1, "popTime")

    expect_s3_class(out2, "data.frame")
    expect_s3_class(out2, "data.table")
    expect_s3_class(out2, "popTime")
    expect_equal(attr(out2, "exposure"), "z")


    expect_s3_class(out3, "data.frame")
    expect_s3_class(out3, "data.table")
    expect_s3_class(out3, "popTime")

    expect_s3_class(out4, "data.frame")
    expect_s3_class(out4, "data.table")
    expect_s3_class(out4, "popTime")
    expect_equal(attr(out4, "exposure"), "z")

})


test_that("plot methods-no error in popTime with data.frame or data.table", {
    p1 <- try(plot(out1,
                   add.case.series = F))
    p2 <- try(plot(out1,
                   add.case.series = T))
    p3 <- try(plot(out1,
                   add.case.series = T,
                   add.base.series = T,
                   ratio = 1,
                   comprisk = T))

    # all three types plotted
    p5 <- try(plot(out1,
                   add.case.series = T,
                   add.base.series = T,
                   add.competing.event = T,
                   ratio = 1,
                   comprisk = T,
                   legend = T,
                   theme.params = list(legend.position = "top")))

    # change points and labels
    p6 <- try(plot(out1,
                   add.case.series = TRUE,
                   add.base.series = TRUE,
                   add.competing.event = TRUE,
                   ratio = 1,
                   comprisk = TRUE,
                   legend = TRUE,
                   case.params = list(mapping = aes(x = time, y = yc,
                                                    color = "Relapse",
                                                    fill = "Relapse")),
                   base.params = list(mapping = aes(x = time, y = ycoord,
                                                    color = "Base series",
                                                    fill = "Base series")),
                   competing.params = list(mapping = aes(x = time, y = yc,
                                                         color = "Competing event",
                                                         fill = "Competing event")),
                   fill.params = list(name = "Legend Name",
                                        breaks = c("Relapse", "Base series",
                                                   "Competing event"),
                                        values = c("Relapse" = "blue",
                                                   "Competing event" = "yellow",
                                                   "Base series" = "orange")),
                   color.params = list(name = "Legend Name",
                                      breaks = c("Relapse", "Base series",
                                                 "Competing event"),
                                      values = c("Relapse" = "blue",
                                                 "Competing event" = "yellow",
                                                 "Base series" = "orange")),
                   theme.params = list(legend.position = "right")))


    # change points and labels and exposure stratified
    p7 <- try(plot(out2,
                   add.case.series = TRUE,
                   add.base.series = TRUE,
                   add.competing.event = TRUE,
                   ratio = 1,
                   comprisk = TRUE,
                   legend = TRUE,
                   case.params = list(mapping = aes(x = time, y = yc,
                                                    color = "Relapse",
                                                    fill = "Relapse")),
                   base.params = list(mapping = aes(x = time, y = ycoord,
                                                    color = "Base series",
                                                    fill = "Base series")),
                   competing.params = list(mapping = aes(x = time, y = yc,
                                                         color = "Competing event",
                                                         fill = "Competing event")),
                   fill.params = list(name = "Legend Name",
                                        breaks = c("Relapse", "Base series",
                                                   "Competing event"),
                                        values = c("Relapse" = "blue",
                                                   "Competing event" = "yellow",
                                                   "Base series" = "orange")),
                   color.params = list(name = "Legend Name",
                                      breaks = c("Relapse", "Base series",
                                                 "Competing event"),
                                      values = c("Relapse" = "blue",
                                                 "Competing event" = "yellow",
                                                 "Base series" = "orange")),
                   casebase.theme = FALSE))

    expect_false(inherits(p1, "try-error"))
    expect_false(inherits(p2, "try-error"))
    expect_false(inherits(p3, "try-error"))
    expect_false(inherits(p5, "try-error"))
    expect_false(inherits(p6, "try-error"))
    expect_false(inherits(p7, "try-error"))

})
