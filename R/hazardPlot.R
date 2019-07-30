
pacman::p_load(reprex)
# simulated non-PH example ------------------------------------------------

pacman::p_load(simsurv) # to simulate data
pacman::p_load(rstpm2) # package that can fit non-PH model
pacman::p_load(casebase)
pacman::p_version(casebase)
pacman::p_load(splines) # for bsplines


covs <- data.frame(id = 1:5000, trt = rbinom(5000, 1, 0.5))
simdat <- simsurv(dist = "weibull", lambdas = 0.1, gammas = 1.5, betas = c(trt = -0.5),
                  x = covs, tde = c(trt = 0.15), tdefunction = "log", maxt = 5)
simdat <- merge(simdat, covs)
head(simdat)

mod_tvc <- rstpm2::stpm2(Surv(eventtime, status) ~ trt, 
                         data = simdat, tvc = list(trt = 1))

mod_cb <- casebase::fitSmoothHazard(status ~ trt + nsx(log(eventtime), df = 3) + 
                                      trt:nsx(log(eventtime),df=1),
                                    time = "eventtime",
                                    data = simdat, 
                                    ratio = 100)
summary(mod_tvc)
summary(mod_cb)




# plot hazard function-rstpm2 ----------------------------------------------------

plot(mod_tvc, newdata = data.frame(trt = 1), type = "hr", 
     var = "trt", ylim = c(0,1), ci = TRUE, rug = FALSE,
     main = "Time dependent hazard ratio",
     ylab = "Hazard ratio", xlab = "Time")
mod_ph <- rstpm2::stpm2(Surv(eventtime, status) ~ trt, 
                        data = simdat)
plot(mod_ph,  newdata = data.frame(trt = 0), type = "hr", 
     var = "trt", ylim = c(0,1), add = TRUE, ci = FALSE, lty = 2)

rstpm2::plot.aft.base()
pacman::p_functions("rstpm2")
rstpm2:::plot.stpm2.base
rstpm2:::predict.stpm2.base



# plot-hazard function casebase -------------------------------------------

pacman::p_load(simsurv) # to simulate data
pacman::p_load(casebase)
pacman::p_version(casebase)
pacman::p_load(splines) # for bsplines


covs <- data.frame(id = 1:5000, trt = rbinom(5000, 1, 0.5))
simdat <- simsurv(dist = "weibull", lambdas = 0.1, gammas = 1.5, betas = c(trt = -0.5),
                  x = covs, tde = c(trt = 0.15), tdefunction = "log", maxt = 5)
simdat <- merge(simdat, covs)

mod_cb <- casebase::fitSmoothHazard(status ~ trt + ns(log(eventtime), df = 3) + 
                                      trt:ns(log(eventtime),df=1),
                                    time = "eventtime",
                                    data = simdat, 
                                    ratio = 100)
newdat <- data.frame(trt = 1, 
                     eventtime = seq(min(simdat$eventtime),max(simdat$eventtime),0.01),
                     offset = 0)

preds <- predict(mod_cb, newdata = newdat, se.fit = TRUE)

newdat$haz <- preds$fit
newdat$se <- preds$se.fit
newdat$lower <- newdat$haz + qnorm(0.025) * newdat$se
newdat$upper <- newdat$haz + qnorm(0.975) * newdat$se

matplot(newdat$eventtime, exp(newdat$haz), type = "n", xlab = "Time", ylab = "Hazard function")
polygon(c(newdat$eventtime, rev(newdat$eventtime)), 
        c(exp(newdat$lower), rev(exp(newdat$upper))), 
        col = "grey", 
        border = "grey")
lines(newdat$eventtime, exp(newdat$haz), col = 1, lty = 1)




# Make a function to plot hazard function ---------------------------------


function (x, y, newdata = NULL, type = "surv", xlab = NULL, 
          ylab = NULL, line.col = 1, ci.col = "grey", lty = par("lty"), 
          log = "", add = FALSE, ci = !add, rug = !add, var = NULL, 
          exposed = incrVar(var), times = NULL, type.relsurv = c("excess", 
                                                                 "total", "other"), ratetable = survival::survexp.us, 
          rmap, scale = 365.24, ...) 
{
  if (type %in% c("meansurv", "meansurvdiff", "af", "meanhaz", 
                  "meanhr")) {
    return(plot.meansurv(x, times = times, newdata = newdata, 
                         type = type, xlab = xlab, ylab = ylab, line.col = line.col, 
                         ci.col = ci.col, lty = lty, add = add, ci = ci, 
                         rug = rug, exposed = exposed, ...))
  }
  if (is.null(newdata)) 
    stop("newdata argument needs to be specified")
  y <- predict(x, newdata, type = switch(type, fail = "surv", 
                                         margfail = "margsurv", type), var = var, exposed = exposed, 
               grid = !(x@timeVar %in% names(newdata)), se.fit = ci, 
               keep.attributes = TRUE, type.relsurv = type.relsurv, 
               ratetable = ratetable, rmap = rmap, scale = scale)
  if (type %in% c("fail", "margfail")) {
    if (ci) {
      y$Estimate <- 1 - y$Estimate
      lower <- y$lower
      y$lower = 1 - y$upper
      y$upper = 1 - lower
    }
    else y <- structure(1 - y, newdata = attr(y, "newdata"))
  }
  if (is.null(xlab)) 
    xlab <- deparse(x@timeExpr)
  if (is.null(ylab)) 
    ylab <- switch(type, hr = "Hazard ratio", hazard = "Hazard", 
                   surv = "Survival", density = "Density", sdiff = "Survival difference", 
                   hdiff = "Hazard difference", cumhaz = "Cumulative hazard", 
                   loghazard = "log(hazard)", link = "Linear predictor", 
                   meansurv = "Mean survival", meansurvdiff = "Difference in mean survival", 
                   meanhr = "Mean hazard ratio", odds = "Odds", or = "Odds ratio", 
                   margsurv = "Marginal survival", marghaz = "Marginal hazard", 
                   marghr = "Marginal hazard ratio", haz = "Hazard", 
                   fail = "Failure", meanhaz = "Mean hazard", margfail = "Marginal failure", 
                   af = "Attributable fraction", meanmargsurv = "Mean marginal survival", 
                   uncured = "Uncured distribution")
  xx <- attr(y, "newdata")
  xx <- eval(x@timeExpr, xx)
  if (!add) 
    matplot(xx, y, type = "n", xlab = xlab, ylab = ylab, 
            log = log, ...)
  if (ci) {
    polygon(c(xx, rev(xx)), c(y[, 2], rev(y[, 3])), col = ci.col, 
            border = ci.col)
    lines(xx, y[, 1], col = line.col, lty = lty, ...)
  }
  else lines(xx, y, col = line.col, lty = lty, ...)
  if (rug) {
    Y <- x@y
    eventTimes <- Y[Y[, ncol(Y)] == 1, ncol(Y) - 1]
    rug(eventTimes, col = line.col)
  }
  return(invisible(y))
}





















rstpm2:::plot.stpm2.base



pacman::p_load(dplyr)

newdat0 <- newdat %>% filter(trt==0)
newdat1 <- newdat %>% filter(trt==1)

plot(newdat0$eventtime, newdat1$haz / newdat0$haz, type = "l")
plot(mod_tvc,  newdata = data.frame(trt = 0), type = "hr", 
     var = "trt", 
     # ylim = c(0,1), 
     add = TRUE, ci = FALSE, lty = 2)


newdat2 <- data.frame(trt = 1, 
                      eventtime = seq(min(simdat$eventtime),max(simdat$eventtime),0.01),
                      offset = 0)
newdat2$haz <- exp(predict(mod_cb, newdata = newdat2))
plot(newdat2$eventtime, newdat2$haz, type = "p", pch = 19, col = "red")

library(magrittr)
newdat %>% head
rstpm2:::plot.stpm2.base

plot(mod_tvc, newdata = data.frame(trt = 0), type = "hazard", 
     var = "trt", 
     # ylim = c(0,1), 
     ci = FALSE, 
     rug = FALSE,
     # main = "Time dependent hazard ratio",
     # ylab = "Hazard ratio", 
     xlab = "Time")
lines(ttime, exp(haz), type = "l")
