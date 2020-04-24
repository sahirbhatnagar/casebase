# devtools::load_all()
# data("simdat")
# library(splines)
# mod_cb <- casebase::fitSmoothHazard(status ~ trt + ns(log(eventtime), df = 3) +
#                                         trt:ns(log(eventtime),df=1),
#                                     time = "eventtime",
#                                     data = simdat,
#                                     ratio = 100,
#                                     family = "glm")
# par(mfrow = c(1,2))
# results00 <- hazardPlotNew(object = mod_cb, newdata = data.frame(trt = 0),
#                           ci.lvl = 0.95, ci = TRUE, lty = 1, line.col = 1, lwd = 2)
#
# results0 <- hazardPlot(object = mod_cb, newdata = data.frame(trt = 0),
#            ci.lvl = 0.95, ci = TRUE, lty = 1, line.col = 1, lwd = 2)
#
#
# results00$predictedloghazard/results0$predictedloghazard
#
# head(results0)
# hazardPlot(object = mod_cb, newdata = data.frame(trt = 1), ci = TRUE,
#            ci.lvl = 0.95, add = TRUE, lty = 2, line.col = 2, lwd = 2)
# legend("topleft", c("trt=0","trt=1"),lty=1:2,col=1:2,bty="y", lwd = 2)
#
#
#
# family(mod_cb)
#
# .clean_terms <- function(x) {
#     # get positions of variable names and see if we have
#     # a suffix for certain values
#     cleaned.pos <- regexpr(pattern = "(\\s|\\[)", x)
#
#     # position "-1" means we only had variable name, no suffix
#     replacers <- which(cleaned.pos == -1)
#     # replace -1 with number of chars
#     cleaned.pos[replacers] <- nchar(x)[replacers]
#
#     # get variable names only
#     x <- trimws(substr(x, 0, cleaned.pos))
#
#     # be sure to remove any brackets
#     sub("[", "", x, fixed = TRUE)
# }
#
# .clean_terms(mod_cb$formula)
# mod_cb$formula
#
# plotHR <- function (model , terms = 1 , se = TRUE , col.se = "#DEEBF7",
#                     rug = "density" , col.dens = grey(0.8) ,
#                     ref = FALSE , col.ref = "#000000" , lty.ref = 1 , lwd.ref = 2 ,
#                     xlab = "" , ylab = "Hazard Ratio" , main = NULL ,
#                     xlim = NULL , ylim = NULL,
#                     col.term = "#08519C", lwd.term = 3, digits = 1,
#                     cex = 1 , bty = "n" , axes = TRUE , ...) {
#
#
#
#
#     # set plotting parameters
#     par(las = 1 , cex = cex)
#
#
#     ## extract the names of all model covariates
#     all.labels <- attr(model$terms , "term.labels")
#
#     # remove 'strata()' / 'factor()' / "pspline( , df = a )"
#     all.labels <- sub ( "pspline[(]([a-zA-Z._0-9]*)(|[, ]*(df|nterm)[ ]*=[ ]*[0-9a-zA-Z_]+)[)]" , "\\1" , all.labels)
#     all.labels <- sub ( "factor[(]([a-zA-Z._0-9]*)[)]" , "\\1" , all.labels)
#     all.labels <- sub ( "strata[(]([a-zA-Z._0-9]*)[)]" , "\\1" , all.labels)
#
#     # The cph in the Design package uses rcs instead of pspline
#     all.labels <- sub ( "rcs[(]([a-zA-Z._0-9]*)(|([, ]*([0-9]+|cut[ ]*=[ ]*c[(][0-9a-zA-Z_, ]+[)])))[)]" , "\\1" , all.labels)
#
#     # Allow the term searched for be a string
#     if (is.character(terms)){
#         terms <- grep(terms, all.labels)
#         if(length(terms) == 0){
#             stop(paste("Could not find term:", terms))
#         }
#     }
#
#     # pick the name of the main term which is goint to be plotted
#     term.label <- all.labels[terms]
#     if (is.na(term.label)){
#         stop(paste("Term", terms, "not found"))
#     }
#
#     # Remove interaction terms since the data can't be found, ex. male_gender:prosthesis
#     terms_with_interaction <- grep("[_.a-zA-Z0-9]+:[_.a-zA-Z0-9]+", all.labels)
#     if(length(terms_with_interaction)>0){
#         all.labels <- all.labels[!(terms_with_interaction == 1:length(all.labels))]
#     }
#
#     ## extract data from model;
#     # only the covariates really used in the model
#     # only complete covariate records (only used in the model anyway)
#     # 'as.data.frame()' and 'names()' have to be explicitly specified in case of a univariate model
#     data <- eval(model$call$data)
#     data <- as.data.frame(na.exclude(data[ , all.labels]))
#     names(data) <- all.labels
#
#
#     ## get the quartiles of the main term
#     quantiles <- quantile(data[,term.label] , probs = c(0.025,0.25,0.50,0.75,0.975))
#
#
#     ### _______________ the smooth term prediction ____________
#     ## calculate the HR for all the covariate values found in the dataset
#     if(length(grep("cph", model)) > 0){
#         # If this is a cph model then don't exclude the na values
#         term <- predict (model , type="terms" , se.fit = TRUE , terms = terms, expand.na=FALSE, na.action=na.pass)
#     }else{
#         term <- predict (model , type="terms" , se.fit = TRUE , terms = terms)
#     }
#
#     # attach the smooth fit for the HR ('fit') and the CI's to the dataset
#     data$fit <- term$fit
#     data$ucl <- term$fit + 1.96 * term$se.fit
#     data$lcl <- term$fit - 1.96 * term$se.fit
#     ### _______________ this is the main extraction __________
#
#
#
#     ### _______________ what now follows is the graphical manipulation of the plots ______
#
#     # plot empty plot with coordinatesystem and labels
#     plot( data$fit ~ data[ , term.label] , xlim = xlim , ylim = ylim , xlab = xlab , ylab = ylab , main = main, type = "n" , axes = FALSE )
#
#
#     # sorting indices for ploting
#     i.backw <- order(data[,term.label] , decreasing = TRUE)
#     i.forw <- order(data[,term.label])
#
#
#     # plot CI as polygon shade - if 'se = TRUE' (default)
#     if (se) {
#         x.poly <- c(data[,term.label][i.forw] , data[,term.label][i.backw])
#         y.poly <- c(data$ucl[i.forw] , data$lcl[i.backw])
#         polygon(x.poly , y.poly , col = col.se , border = NA)
#     }
#
#
#     ### _______________ rug = "density" ____________________________________________________
#     ### density plot at bottom of panel if rug = "density" in function call
#
#     if (rug == "density") {
#
#         # calculate the coordinates of the density function
#
#         density <- density( data[,term.label] )
#
#         # the height of the densityity curve
#
#         max.density <- max(density$y)
#
#         # Get the boundaries of the plot to
#         # put the density polygon at the x-line
#
#         plot_coordinates <- par("usr")
#
#         # get the "length" and range of the y-axis
#         y.scale <- plot_coordinates[4] - plot_coordinates[3]
#
#         # transform the y-coordinates of the density
#         # to the lower 10% of the plotting panel
#
#         density$y <- (0.1 * y.scale / max.density) * density$y + plot_coordinates[3]
#
#         ## plot the polygon
#
#         polygon( density$x , density$y , border = F , col = col.dens) }
#
#
#     # plot white lines (background color) for 2.5%tile, 1Q, median, 3Q and 97.5%tile through confidence shade and density plot
#     #  axis( side = 1 , at = quantiles , labels = FALSE , lwd = 0 , col.ticks = "white"  , lwd.ticks = 2 , tck = 1 )
#     # plot white line ( background color for HR = 1 reference
#     if (ref) {
#         axis( side = 2 , at = 0 , labels = FALSE , lwd = 0 , col.ticks = col.ref  , lwd.ticks = lwd.ref , lty.ticks = lty.ref , tck = 1 )}
#
#
#     ### _______________ rug = "ticks" ____________________________________________________
#     ### rug plot if "ticks" is specified in function call
#     if (rug == "ticks") {
#         # rugs at datapoints
#         axis(side = 1 , line = -1.2 , at = jitter(data[,term.label]) , labels = F , tick = T , tcl = 0.8 , lwd.ticks = 0.1 , lwd = 0)
#         # rugs and labels at 1Q, median and 3Q
#         axis(side = 1 , line = -1.0 , at = fivenum(data[,term.label])[2:4], lwd = 0 , tick = T, tcl = 1.2 , lwd.ticks = 1 , col.ticks = "black" , labels = c("Quartile 1","Median","Quartile 3"), cex.axis = 0.7, col.axis = "black" , padj = -2.8)
#         axis(side = 1 , line = 0.0 , at = fivenum(data[,term.label])[2:4], lwd = 0 , tick = T, tcl = 0.2 , lwd.ticks = 1 , col.ticks = "black", labels = FALSE)
#     }
#
#
#     # last but not least the main plotting line for the smooth estimate:
#     # ___________ main plot _____________
#     lines(data[,term.label][i.forw] , data$fit[i.forw], col = col.term, lwd = lwd.term)
#     # ___________ main plot _____________
#
#
#     # plot the axes
#     if (axes){
#         axis(side = 1)
#         axis(side = 2 , at = axTicks(2) , label = round( exp(axTicks(2)) , digits = digits))
#     }
#
#
#     # plot a box around plotting panel if specified - not plotted by default
#     box(bty = bty)
# }
#
# data("simdat")
# library(splines)
# library(survival)
# mod_cb <- casebase::fitSmoothHazard(status ~ trt + ns(log(eventtime), df = 3) +
#                                         trt:ns(log(eventtime),df=1),
#                                     time = "eventtime",
#                                     data = simdat,
#                                     ratio = 100,
#                                     family = "glm")
# head(simdat)
# data("veteran")
# head(veteran)
# model <- coxph(Surv(time, status)~trt + karno, data=veteran)
# summary(model)
# dev.off()
# plotHR(model, terms=1, se = "TRUE" , rug = "density" , xlab = "" , ylab = "" ,
#        xlim=c(), ylim=c(), col.term = "#08519C",
#        lwd.term = "2",
#        col.se = "#DEEBF7",
#        cex = "1.5" ,
#        bty = "L" , axes = TRUE)
#
#
# insight::find_predictors(mod_cb)
#
# pacman::p_load_gh("jacob-long/interactions")
# fiti <- lm(mpg ~ hp * wt, data = mtcars)
# head(mtcars)
# sim_slopes(fiti, pred = hp, modx = wt, jnplot = TRUE)
# sim_slopes(mod_cb, pred = hp, modx = wt, jnplot = TRUE)
#
# interact_plot(fiti, pred = hp, modx = wt, interval = TRUE)
#
# head(simdat)
# devtools::load_all("~/git_repositories/casebase/")
#
#
# library(survival)
# data("jasa")
# head(simdat)
# str(veteran)
# Surv(time, status)~trt + karno, data=veteran
# library(splines)
# mod_cb <- casebase::fitSmoothHazard(status ~ bs(eventtime)*factor(trt),
#                                     time = "eventtime",
#                                     data = simdat,
#                                     ratio = 100,
#                                     family = "glm")
# summary(mod_cb)
# # interact_plot(mod_cb, pred = eventtime, modx = trt, interval = TRUE, data = mod_cb$data)
#
# pacman::p_load_gh("leeper/margins")
# mod_cb$originalData %>% head
# mod_cb$originalData$offset <- 0
# tt <- prediction::prediction(mod_cb, data = mod_cb$originalData)
# head(tt)
# tt %>% dim
# ff <- dydx(data = mod_cb$originalData, model = mod_cb, variable = "trt")
# plot(mod_cb$originalData$eventtime, ff$dydx_trt)
# tt$fitted
# tt$se.fitted
# plot(tt$time, tt$fitted)
# margins::marginal_effects(data = mod_cb$originalData, model = mod_cb)  %>% head
#
# margins::margins(data = mod_cb$originalData, model = mod_cb)
#
#
#
#
#
#
# data("birthwt", package="MASS")
# fit <- glm(low ~ age + race + smoke + lwt, data=birthwt, family="binomial")
# pacman::p_load(visreg)
#
# dev.off()
# visreg(mod_cb,
#        xvar = "eventtime",
#        by = "trt",
#        xlab="Mother's weight", ylab="Log odds (low birthweight)")
#
# visreg(fit, "lwt", xlab="Mother's weight", ylab="Log odds (low birthweight)")
# visreg:::setupF
# visreg(fit, "lwt", scale="response", rug=2, xlab="Mother's weight",
#        ylab="P(low birthweight)")
#
#
#
#
#
#
#
#
#
# pacman::p_load(rstpm2) # for brcancer data
# pacman::p_load(magrittr)
# pacman::p_load(visreg)
# devtools::load_all("~/git_repositories/casebase/")
# data("brcancer")
# brcancer <- transform(brcancer, recyear=rectime / 365.24)
# brcancer <- transform(brcancer, hormon=factor(hormon))
# head(brcancer)
# mod_cb <- fitSmoothHazard(censrec ~ hormon*bs(recyear, df = 5),
#                           data = brcancer,
#                           time = "recyear")
# summary(mod_cb)
# # mod_cb$data <- mod_cb$originalData
# # resid(mod_cb)
# # length(resid(mod_cb))
# # stats:::residuals.glm
# # mod_cb$residuals
# # mod_cb$data$offset <- 0
# # mod_cb$originalData$offset <- rnorm(nrow(mod_cb$originalData))
#
#
# par(mfrow=c(1,2))
# tt <- visreg(mod_cb,
#        xvar = "recyear",
#        # type = "contrast",
#        by = "hormon",
#        scale = "linear",
#        cond = list(offset = 0),
#        trans = exp,
#        # data = mod_cb$originalData,
#        plot = T,
#        alpha = 1,
#        overlay = TRUE)
#
# tt <- visreg(mod_cb,
#              xvar = "recyear",
#              type = "contrast",
#              by = "hormon",
#              scale = "linear",
#              cond = list(offset = 0),
#              trans = exp,
#              # data = mod_cb$originalData,
#              plot = T,
#              alpha = 1,
#              overlay = TRUE)
#
# traceback()
# mod_cb$data %>% dim
# mod_cb$originalData %>% dim
# visreg:::setupV
# tt$fit
#
# visreg:::visreg()
#
#
# visreg::visreg
# visreg:::setupV
#
#
#
# # comparing type=contrast in visreg ---------------------------------------
#
# pacman::p_load(rstpm2) # for brcancer data
# pacman::p_load(magrittr)
# pacman::p_load(visreg)
# devtools::load_all("~/git_repositories/casebase/")
# data("brcancer")
# brcancer <- transform(brcancer, recyear=rectime / 365.24)
# brcancer <- transform(brcancer, hormon=factor(hormon))
# mod_cb <- fitSmoothHazard(censrec ~ hormon*bs(recyear, df = 5),
#                           data = brcancer,
#                           time = "recyear")
# summary(mod_cb)
# par(mfrow=c(1,2))
# t0_hazardplot <- hazardPlotNew(object = mod_cb, newdata = data.frame(hormon = "0"), ci = F,
#                             ci.lvl = 0.95, lty = 2, line.col = 2, lwd = 2)
#
# t1_hazardplot <- hazardPlotNew(object = mod_cb, newdata = data.frame(hormon = "1"),
#                                ci.lvl = 0.95, ci = F, lty = 1, line.col = 1, lwd = 2, add = T)
#
# # this uses the augmented casebase dataset which has a non-zero offset.
# # need to set it to 0 here
# tt <- visreg(mod_cb,
#              xvar = "recyear",
#              # type = "contrast",
#              by = "hormon",
#              scale = "linear",
#              cond = list(offset = 0),
#              trans = exp,
#              # data = mod_cb$originalData,
#              plot = T,
#              alpha = 1,
#              overlay = TRUE)
#
# # get the predicted hazards by trt status and order the
# # times so you can plot them
# t1 <- tt$fit[which(tt$fit$hormon==1),]
# t0 <- tt$fit[which(tt$fit$hormon==0),]
# t1_visreg <- t1[order(t1$recyear),]
# t0_visreg <- t0[order(t1$recyear),]
#
#
# # dev.off()
# par(mfrow = c(1,2))
# plot(t0_hazardplot$recyear, t1_hazardplot$predictedhazard / t0_hazardplot$predictedhazard,
#      type = "l", col = 1, lty = 1, ylab = "hazard ratio for hormon=1 vs. hormon=0", lwd = 2, ylim = c(0,1.1))
# abline(a=1, b=0, col = "grey80")
# lines(t0_visreg$recyear, t1_visreg$visregFit / t0_visreg$visregFit,
#       type = "l", col = 2, lty = 2, lwd = 2)
# legend("bottomleft", c("hazardPlot","visreg"),
#        lty=1:2,col=1:2,bty="y", lwd = 2)
#
# # this is type contrast: E(y|x=1, z) - E(y|x=0, z)
# # where z are the other covariates not of interest
# # I thought it would give me just one line, but its giving me two
# kk <- visreg(mod_cb,
#              xvar = "recyear",
#              type = "contrast",
#              by = "hormon",
#              scale = "linear",
#              cond = list(offset = 0),
#              trans = exp,
#              plot = T,
#              alpha = 1,
#              overlay = T)
#
# # visreg:::plot.visreg
# # visreg:::visregOverlayPlot
#
# k1 <- kk$fit[which(kk$fit$hormon==1),]
# k1_visreg <- k1[order(k1$recyear),]
# plot(k1_visreg$recyear, k1_visreg$visregFit,
#       type = "l", col = 2, lty = 2, lwd = 2)
#
#
#
#
#
#
#
# tt <- visreg(mod_cb,
#              xvar = "recyear",
#              type = "contrast",
#              # by = "hormon",
#              scale = "linear", # this is plotted on the scale of the linear predictor.. so its b0 + b1*x + ... which is equal to log(hazard) in casebase
#              cond = list(offset = 0, hormon = 1),
#              trans = exp,
#              # data = mod_cb$originalData,
#              plot = T,
#              alpha = 1,
#              overlay = TRUE)
#
#
#
#
#
# # Setup
# library(survival)
# data(ovarian)
# ovarian$rx <- factor(ovarian$rx)
# head(ovarian)
# str(ovarian)
#
# str(ovarian)
# # Basic
# fit <- coxph(Surv(futime,fustat) ~ age + rx, data=ovarian)
# summary(fit)
# visreg(fit, "age")
# visreg(fit, "age", type="contrast")
# visreg(fit, "rx")
# visreg(fit, "age", trans=exp)
# visreg(fit, "age", trans=exp, ylim=c(0,20))
#
# # Interaction
# fit <- coxph(Surv(futime,fustat)~age*rx,data=ovarian)
# par(mfrow=c(1,2))
# visreg(fit, "age", cond=list(rx="1"))
# visreg(fit, "age", cond=list(rx="2"))
# visreg(fit,"age",by="rx")
# visreg2d(fit, x="age", y="rx")
#
# # Splines
# require(splines)
# fit <- coxph(Surv(futime,fustat)~ns(age,4)*rx,data=ovarian)
# par(mfrow=c(1,1))
# visreg(fit, "age")
# visreg(fit, xvar = "age", by = "rx", overlay = T, alpha = 1, type = "contrast")
# visreg(fit, "age", type="contrast")
# visreg(fit, "rx", type="contrast")
#
# # Strata
# ovarian$Group <- factor(ovarian$rx)
# fit <- coxph(Surv(futime,fustat) ~ age+strata(Group), data=ovarian)
# visreg(fit, "age")
# visreg(fit, "age", type='contrast')
#
# # Logical
# ovarian$rx <- ovarian$rx == 2
# fit <- coxph(Surv(futime,fustat) ~ age + rx, data=ovarian)
# visreg(fit, "age")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
