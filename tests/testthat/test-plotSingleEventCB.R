context("plot.singleEventCB")
# Uncomment next line to skip tests in non-interactive session
skip_if_not_installed("glmnet")
skip_if_not_installed("gbm")
skip_if_not_installed("splines")
skip_if_not_installed("visreg")

library(splines)
library(visreg)
data("simdat")
data("brcancer")

mod_glm <- casebase::fitSmoothHazard(status ~ trt + ns(log(eventtime), df = 3) +
                                         trt:ns(log(eventtime), df = 1),
                                     time = "eventtime",
                                     data = simdat[sample(nrow(simdat),
                                                          size = 200), ],
                                     ratio = 1,
                                     family = "glm")

mod_brcancer <- fitSmoothHazard(cens ~ ns(time, df = 3) * tgrade,
                                data = brcancer,
                                time = "time")

test_that("no error in plot method for singleEventCB objects - hazard function", {
    outglm <- try(plot(mod_glm,
                       hazard.params = list(xvar = "eventtime",
                                            by = "trt",
                                            alpha = 0.05,
                                            ylab = "Hazard")),
                  silent = TRUE)

    expect_false(inherits(outglm, "try-error"))
})


test_that("no error in plot method for singleEventCB objects - hazard ratio with and without confidence interval", {

    newtime <- quantile(mod_glm[["originalData"]][[mod_glm[["timeVar"]]]],
                        probs = seq(0.01, 0.99, 0.01))

    # reference category
    newdata <- data.frame(trt = 0, eventtime = newtime)


    # check for correctly specified var and exposed args
    expect_error(plot(mod_glm,
                      type = "hr",
                      newdata = newdata,
                      increment = 1,
                      xvar = "eventtime",
                      ci = TRUE,
                      rug = TRUE))

    # var not found in newdata
    expect_error(plot(mod_glm,
                      type = "hr",
                      newdata = newdata,
                      var = "treatment",
                      increment = 1,
                      xvar = "eventtime",
                      ci = TRUE,
                      rug = TRUE))

    # exposed takes priority
    expect_warning(plot(mod_glm,
                      type = "hr",
                      newdata = newdata,
                      exposed = function(data) transform(data,
                                                         trt = 1),
                      var = "treatment",
                      increment = 1,
                      xvar = "eventtime",
                      ci = TRUE,
                      rug = TRUE))

    # exposed isnt a function and var not found in newdata
    expect_error(plot(mod_glm,
                        type = "hr",
                        newdata = newdata,
                        exposed = transform(newdata, trt = 1),
                        var = "treatment",
                        increment = 1,
                        xvar = "eventtime",
                        ci = TRUE,
                        rug = TRUE))

    # exposed isnt a function. use var in newdata
    expect_warning(plot(mod_glm,
                        type = "hr",
                        newdata = newdata,
                        exposed = transform(newdata, trt = 1),
                        var = "trt",
                        increment = 1,
                        xvar = "eventtime",
                        ci = TRUE,
                        rug = TRUE))

    # check that supplying only one row to newdata works
    # and that the dataset is returned invisibly
    expect_invisible(plot(mod_brcancer,
                          type = "hr",
                          newdata = brcancer[1, ],
                          var = "tgrade",
                          increment = 1,
                          xvar = "time",
                          ci = T,
                          rug = T))

    # check that supplying only one row to newdata works
    # and that the dataset is returned invisibly
    expect_invisible(plot(mod_brcancer,
                          type = "hr",
                          newdata = brcancer[1, ],
                          var = "tgrade",
                          increment = 1,
                          xvar = "time",
                          ci = F,
                          rug = T, cex = 5))

    # check that supplying only one row to newdata works
    # and that the dataset is returned invisibly
    # when CI=F and other plot params are passed to ...
    expect_invisible(plot(mod_brcancer,
                          type = "hr",
                          newdata = brcancer[1, ],
                          var = "tgrade",
                          increment = 1,
                          xvar = "time",
                          ci = F,
                          rug = T, cex = 5, pch = 22)
    )


    outglm_hr <- try(plot(mod_glm,
                          type = "hr",
                          newdata = newdata,
                          var = "trt",
                          increment = 1,
                          xvar = "eventtime",
                          ci = TRUE,
                          rug = TRUE),
                     silent = TRUE)

    outglm_hr_noci <- try(plot(mod_glm,
                               type = "hr",
                               newdata = newdata,
                               var = "trt",
                               increment = 1,
                               xvar = "eventtime",
                               ci = FALSE,
                               rug = TRUE),
                          silent = TRUE)

    # using the exposed argument instead
    outglm_hr_exposed <- try(plot(mod_glm,
                                  type = "hr",
                                  newdata = newdata,
                                  exposed = function(data) transform(data,
                                                                     trt = 1),
                                  xvar = "eventtime",
                                  ci = TRUE,
                                  rug = TRUE),
                             silent = TRUE)

    expect_false(inherits(outglm_hr, "try-error"))
    expect_false(inherits(outglm_hr_noci, "try-error"))
    expect_false(inherits(outglm_hr_exposed, "try-error"))
})
