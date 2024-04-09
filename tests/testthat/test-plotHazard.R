context("plotHazard")
set.seed(12345)

# CRAN skip atlas check fix
testthat::skip_if(grepl(pattern = "atlas", sessionInfo()$BLAS,
                        ignore.case = TRUE))

skip_if_not_installed("glmnet")
skip_if_not_installed("mgcv")
skip_if_not_installed("gbm")
skip_if_not_installed("splines")

library(splines)
data("simdat")
mod_glm <- casebase::fitSmoothHazard(status ~ trt + ns(log(eventtime), df = 3) +
                                        trt:ns(log(eventtime), df = 1),
                                    time = "eventtime",
                                    data = simdat,
                                    ratio = 10,
                                    family = "glm")

mod_gam <- casebase::fitSmoothHazard(status ~ trt + s(log(eventtime)) +
                                         s(trt, log(eventtime)),
                                     time = "eventtime",
                                     data = simdat,
                                     ratio = 10,
                                     family = "gam")

mod_glmnet <- casebase::fitSmoothHazard(status ~ trt + ns(log(eventtime),
                                                          df = 3) +
                                        trt:ns(log(eventtime), df = 1),
                                    time = "eventtime",
                                    data = simdat,
                                    ratio = 10,
                                    family = "glmnet")

# mod_gbm <- casebase::fitSmoothHazard(status ~ trt + eventtime,
#                                      time = "eventtime",
#                                      interaction.depth = 2,
#                                      data = simdat,
#                                      ratio = 10,
#                                      family = "gbm")


test_that("no error in hazardPlot for glm, gam, glmnet using default times", {
    # outgbm <- try(hazardPlot(object = mod_gbm, newdata = data.frame(trt = 0),
    #                          add = FALSE, ci.lvl = 0.95, ci = FALSE, lty = 1,
    #                          line.col = 1, lwd = 2))
    outglm <- try(hazardPlot(object = mod_glm, newdata = data.frame(trt = 0),
                             add = FALSE, ci.lvl = 0.95, ci = FALSE, lty = 1,
                             line.col = 2, lwd = 2))
    outgam <- try(hazardPlot(object = mod_gam, newdata = data.frame(trt = 0),
                             add = FALSE, ci.lvl = 0.95, ci = FALSE, lty = 1,
                             line.col = 3, lwd = 2))
    outglmnet <- try(hazardPlot(object = mod_glmnet,
                                newdata = data.frame(trt = 0),
                                add = FALSE, s = "lambda.min", ci.lvl = 0.95,
                                ci = FALSE, lty = 1, line.col = 4, lwd = 2))

    # expect_false(inherits(outgbm, "try-error"))
    expect_false(inherits(outglm, "try-error"))
    expect_false(inherits(outgam, "try-error"))
    expect_false(inherits(outglmnet, "try-error"))
})

test_that("no error in hazardPlot for glm, gam, glmnet using user-defined times", {
    # outgbm <- try(hazardPlot(object = mod_gbm, newdata = data.frame(trt = 0),
    #                          add = FALSE, times = runif(10, 0, 3),
    #                          ci.lvl = 0.95, ci = FALSE, lty = 1, line.col = 1,
    #                          lwd = 2))
    outglm <- try(hazardPlot(object = mod_glm, newdata = data.frame(trt = 0),
                             add = FALSE, times = runif(10, 0, 3),
                             ci.lvl = 0.95, ci = TRUE, lty = 1, line.col = 2,
                             lwd = 2))
    outgam <- try(hazardPlot(object = mod_gam, newdata = data.frame(trt = 0),
                             add = FALSE, times = runif(10, 0, 3),
                             ci.lvl = 0.95, ci = TRUE, lty = 1, line.col = 3,
                             lwd = 2))
    outglmnet <- try(hazardPlot(object = mod_glmnet,
                                newdata = data.frame(trt = 0), add = FALSE,
                                s = "lambda.min", times = runif(10, 0, 3),
                                ci.lvl = 0.95, ci = FALSE, lty = 1,
                                line.col = 4, lwd = 2))

    # expect_false(inherits(outgbm, "try-error"))
    expect_false(inherits(outglm, "try-error"))
    expect_false(inherits(outgam, "try-error"))
    expect_false(inherits(outglmnet, "try-error"))
})
