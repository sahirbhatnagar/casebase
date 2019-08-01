context("plotHazard function")

data("simdat")
mod_glm <- casebase::fitSmoothHazard(status ~ trt + ns(log(eventtime), df = 3) +
                                        trt:ns(log(eventtime),df=1),
                                    time = "eventtime",
                                    data = simdat,
                                    ratio = 10,
                                    family = "glm")

mod_gam <- casebase::fitSmoothHazard(status ~ trt + s(log(eventtime)) +
                                         s(trt,log(eventtime)) ,
                                     time = "eventtime",
                                     data = simdat,
                                     ratio = 10,
                                     family = "gam")

mod_glmnet <- casebase::fitSmoothHazard(status ~ trt + nsx(log(eventtime), df = 3) +
                                        trt:nsx(log(eventtime),df=1),
                                    time = "eventtime",
                                    data = simdat,
                                    ratio = 10,
                                    family = "glmnet")

mod_gbm <- casebase::fitSmoothHazard(status ~ trt + log(eventtime),
                                     time = "eventtime",
                                     interaction.depth = 2,
                                     data = simdat,
                                     ratio = 10,
                                     family = "gbm")
# par(mfrow=c(1,2))
# hazardPlot(object = mod_gbm, newdata = data.frame(trt = 0),add = FALSE,
#            ci.lvl = 0.95, ci = F, lty = 1, line.col = 1, lwd = 2)
# hazardPlot(object = mod_glm, newdata = data.frame(trt = 0),add = T,
#            ci.lvl = 0.95, ci = FALSE, lty = 1, line.col = 2, lwd = 2)
# hazardPlot(object = mod_gam, newdata = data.frame(trt = 0),add = TRUE,
#            ci.lvl = 0.95, ci = F, lty = 1, line.col = 3, lwd = 2)
# hazardPlot(object = mod_glmnet, newdata = data.frame(trt = 0),add = TRUE,s = "lambda.min",
#            ci.lvl = 0.95, ci = F, lty = 1, line.col = 4, lwd = 2)
# legend("topleft", c("GBM","GLM","GAM","GLMNET"),
#        lty=1,col=1:4,bty="y", lwd = 2)
#
#
# hazardPlot(object = mod_gbm, newdata = data.frame(trt = 1),add = FALSE,
#            ci.lvl = 0.95, ci = F, lty = 1, line.col = 1, lwd = 2)
# hazardPlot(object = mod_glm, newdata = data.frame(trt = 1),add = T,
#            ci.lvl = 0.95, ci = FALSE, lty = 1, line.col = 2, lwd = 2)
# hazardPlot(object = mod_gam, newdata = data.frame(trt = 1),add = TRUE,
#            ci.lvl = 0.95, ci = F, lty = 1, line.col = 3, lwd = 2)
# hazardPlot(object = mod_glmnet, newdata = data.frame(trt = 1),add = TRUE,s = "lambda.min",
#            ci.lvl = 0.95, ci = F, lty = 1, line.col = 4, lwd = 2)
# legend("topleft", c("GBM","GLM","GAM","GLMNET"),
#        lty=1,col=1:4,bty="y", lwd = 2)
#
# par(mfrow=c(2,2))
# hazardPlot(object = mod_gbm, newdata = data.frame(trt = 0),add = FALSE,
#            ci.lvl = 0.95, ci = F, lty = 1, line.col = 1, lwd = 2, main = "GBM")
# hazardPlot(object = mod_gbm, newdata = data.frame(trt = 1),add = T,
#            ci.lvl = 0.95, ci = F, lty = 2, line.col = 2, lwd = 2)
# legend("topleft", c("trt=0","trt=1"),lty=1:2,col=1:2,bty="y", lwd = 2)
#
# hazardPlot(object = mod_glm, newdata = data.frame(trt = 0),add = FALSE,
#            ci.lvl = 0.95, ci = F, lty = 1, line.col = 1, lwd = 2, main = "GLM")
# hazardPlot(object = mod_glm, newdata = data.frame(trt = 1),add = T,
#            ci.lvl = 0.95, ci = F, lty = 2, line.col = 2, lwd = 2)
# legend("topleft", c("trt=0","trt=1"),lty=1:2,col=1:2,bty="y", lwd = 2)
#
# hazardPlot(object = mod_gam, newdata = data.frame(trt = 0),add = FALSE,
#            ci.lvl = 0.95, ci = F, lty = 1, line.col = 1, lwd = 2, main = "GAM")
# hazardPlot(object = mod_gam, newdata = data.frame(trt = 1),add = T,
#            ci.lvl = 0.95, ci = F, lty = 2, line.col = 2, lwd = 2)
# legend("topleft", c("trt=0","trt=1"),lty=1:2,col=1:2,bty="y", lwd = 2)
#
# hazardPlot(object = mod_glmnet, newdata = data.frame(trt = 0),add = FALSE,s = "lambda.min",
#            ci.lvl = 0.95, ci = F, lty = 1, line.col = 1, lwd = 2, main = "GLMNET")
# hazardPlot(object = mod_glmnet, newdata = data.frame(trt = 1),add = T,s = "lambda.min",
#            ci.lvl = 0.95, ci = F, lty = 2, line.col = 2, lwd = 2)
# legend("topleft", c("trt=0","trt=1"),lty=1:2,col=1:2,bty="y", lwd = 2)
# dev.off()
#
# hazardPlot(object = mod_cb, newdata = data.frame(trt = 1), ci = TRUE,
#            ci.lvl = 0.95, add = TRUE, lty = 2, line.col = 2, lwd = 2)


test_that("no error in hazardPlot for glm, gbm, gam, glmnet using default times", {
    outgbm <- try(hazardPlot(object = mod_gbm, newdata = data.frame(trt = 0),add = FALSE,
                           ci.lvl = 0.95, ci = F, lty = 1, line.col = 1, lwd = 2))
    outglm <- try(hazardPlot(object = mod_glm, newdata = data.frame(trt = 0),add = F,
                           ci.lvl = 0.95, ci = FALSE, lty = 1, line.col = 2, lwd = 2))
    outgam <- try(hazardPlot(object = mod_gam, newdata = data.frame(trt = 0),add = F,
                           ci.lvl = 0.95, ci = F, lty = 1, line.col = 3, lwd = 2))
    outglmnet <- try(hazardPlot(object = mod_glmnet, newdata = data.frame(trt = 0),add = F,
                                s = "lambda.min",
                           ci.lvl = 0.95, ci = F, lty = 1, line.col = 4, lwd = 2))

    expect_false(inherits(outgbm, "try-error"))
    expect_false(inherits(outglm, "try-error"))
    expect_false(inherits(outgam, "try-error"))
    expect_false(inherits(outglmnet, "try-error"))
})

test_that("no error in hazardPlot for glm, gbm, gam, glmnet using user-defined times", {
    outgbm <- try(hazardPlot(object = mod_gbm, newdata = data.frame(trt = 0),add = FALSE,
                             times = runif(10,0,3),
                             ci.lvl = 0.95, ci = F, lty = 1, line.col = 1, lwd = 2))
    outglm <- try(hazardPlot(object = mod_glm, newdata = data.frame(trt = 0),add = F,
                             times = runif(10,0,3),
                             ci.lvl = 0.95, ci = T, lty = 1, line.col = 2, lwd = 2))
    outgam <- try(hazardPlot(object = mod_gam, newdata = data.frame(trt = 0),add = F,
                             times = runif(10,0,3),
                             ci.lvl = 0.95, ci = T, lty = 1, line.col = 3, lwd = 2))
    outglmnet <- try(hazardPlot(object = mod_glmnet, newdata = data.frame(trt = 0),add = F,
                                s = "lambda.min",
                                times = runif(10,0,3),
                                ci.lvl = 0.95, ci = F, lty = 1, line.col = 4, lwd = 2))

    expect_false(inherits(outgbm, "try-error"))
    expect_false(inherits(outglm, "try-error"))
    expect_false(inherits(outgam, "try-error"))
    expect_false(inherits(outglmnet, "try-error"))
})
