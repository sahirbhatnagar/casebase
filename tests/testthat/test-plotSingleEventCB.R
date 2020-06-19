context("plot.singleEventCB")
# Uncomment next line to skip tests in non-interactive session
# skip_if_not(interactive())
# skip_if_not_installed("glmnet")
# skip_if_not_installed("mgcv")
# skip_if_not_installed("gbm")
skip_if_not_installed("splines")
skip_if_not_installed("visreg")

library(splines)
library(visreg)
data("simdat")
mod_glm <- casebase::fitSmoothHazard(status ~ trt + ns(log(eventtime), df = 3) +
                                         trt:ns(log(eventtime),df=1),
                                     time = "eventtime",
                                     data = simdat[sample(1:nrow(simdat), size = 200),],
                                     ratio = 1,
                                     family = "glm")

test_that("no error in plot method for singleEventCB objects - hazard function", {
    outglm <- try(plot(mod_glm,
                       hazard.params = list(xvar = "eventtime",
                                            by = "trt",
                                            alpha = 0.05,
                                            ylab = "Hazard")))

    expect_false(inherits(outglm, "try-error"))
})


test_that("no error in plot method for singleEventCB objects - hazard ratio", {

    newtime <- quantile(mod_glm[["originalData"]][[mod_glm[["timeVar"]]]], probs = seq(0.01, 0.99, 0.01))

    # reference category
    newdata <- data.frame(trt = 0, eventtime = newtime)

    outglm <- try(plot(mod_glm,
                       hazard.params = list(xvar = "eventtime",
                                            by = "trt",
                                            alpha = 0.05,
                                            ylab = "Hazard")))

    outglm_hr <- try(plot(mod_glm,
                          type = "hr",
                          newdata = newdata,
                          var = "trt",
                          increment = 1,
                          xvar = "eventtime",
                          ci = T,
                          rug = T))

    #using the exposed argument instead
    outglm_hr_exposed <- try(plot(mod_glm,
                                  type = "hr",
                                  newdata = newdata,
                                  exposed = function(data) transform(data, trt = 1),
                                  xvar = "eventtime",
                                  ci = T,
                                  rug = T))

    expect_false(inherits(outglm, "try-error"))
    expect_false(inherits(outglm_hr, "try-error"))
    expect_false(inherits(outglm_hr_exposed, "try-error"))
})

