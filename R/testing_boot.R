plotHR <- function (model , terms = 1 , se = TRUE , col.se = "#DEEBF7",
                    rug = "density" , col.dens = grey(0.8) ,
                    ref = FALSE , col.ref = "#000000" , lty.ref = 1 , lwd.ref = 2 ,
                    xlab = "" , ylab = "Hazard Ratio" , main = NULL ,
                    xlim = NULL , ylim = NULL,
                    col.term = "#08519C", lwd.term = 3, digits = 1,
                    cex = 1 , bty = "n" , axes = TRUE , ...) {




    # set plotting parameters
    par(las = 1 , cex = cex)


    ## extract the names of all model covariates
    all.labels <- attr(model$terms , "term.labels")

    # remove 'strata()' / 'factor()' / "pspline( , df = a )"
    all.labels <- sub ( "pspline[(]([a-zA-Z._0-9]*)(|[, ]*(df|nterm)[ ]*=[ ]*[0-9a-zA-Z_]+)[)]" , "\\1" , all.labels)
    all.labels <- sub ( "factor[(]([a-zA-Z._0-9]*)[)]" , "\\1" , all.labels)
    all.labels <- sub ( "strata[(]([a-zA-Z._0-9]*)[)]" , "\\1" , all.labels)

    # The cph in the Design package uses rcs instead of pspline
    all.labels <- sub ( "rcs[(]([a-zA-Z._0-9]*)(|([, ]*([0-9]+|cut[ ]*=[ ]*c[(][0-9a-zA-Z_, ]+[)])))[)]" , "\\1" , all.labels)

    # Allow the term searched for be a string
    if (is.character(terms)){
        terms <- grep(terms, all.labels)
        if(length(terms) == 0){
            stop(paste("Could not find term:", terms))
        }
    }

    # pick the name of the main term which is goint to be plotted
    term.label <- all.labels[terms]
    if (is.na(term.label)){
        stop(paste("Term", terms, "not found"))
    }

    # Remove interaction terms since the data can't be found, ex. male_gender:prosthesis
    terms_with_interaction <- grep("[_.a-zA-Z0-9]+:[_.a-zA-Z0-9]+", all.labels)
    if(length(terms_with_interaction)>0){
        all.labels <- all.labels[!(terms_with_interaction == 1:length(all.labels))]
    }

    ## extract data from model;
    # only the covariates really used in the model
    # only complete covariate records (only used in the model anyway)
    # 'as.data.frame()' and 'names()' have to be explicitly specified in case of a univariate model
    data <- eval(model$call$data)
    data <- as.data.frame(na.exclude(data[ , all.labels]))
    names(data) <- all.labels


    ## get the quartiles of the main term
    quantiles <- quantile(data[,term.label] , probs = c(0.025,0.25,0.50,0.75,0.975))


    ### _______________ the smooth term prediction ____________
    ## calculate the HR for all the covariate values found in the dataset
    if(length(grep("cph", model)) > 0){
        # If this is a cph model then don't exclude the na values
        term <- predict (model , type="terms" , se.fit = TRUE , terms = terms, expand.na=FALSE, na.action=na.pass)
    }else{
        term <- predict (model , type="terms" , se.fit = TRUE , terms = terms)
    }

    term$fit %>% head
    summary(model)

    head(data)

    browser()

    # attach the smooth fit for the HR ('fit') and the CI's to the dataset
    data$fit <- term$fit
    data$ucl <- term$fit + 1.96 * term$se.fit
    data$lcl <- term$fit - 1.96 * term$se.fit
    ### _______________ this is the main extraction __________



    ### _______________ what now follows is the graphical manipulation of the plots ______

    # plot empty plot with coordinatesystem and labels
    plot( data$fit ~ data[ , term.label] , xlim = xlim , ylim = ylim , xlab = xlab , ylab = ylab , main = main, type = "n" , axes = FALSE )


    # sorting indices for ploting
    i.backw <- order(data[,term.label] , decreasing = TRUE)
    i.forw <- order(data[,term.label])


    # plot CI as polygon shade - if 'se = TRUE' (default)
    if (se) {
        x.poly <- c(data[,term.label][i.forw] , data[,term.label][i.backw])
        y.poly <- c(data$ucl[i.forw] , data$lcl[i.backw])
        polygon(x.poly , y.poly , col = col.se , border = NA)
    }


    ### _______________ rug = "density" ____________________________________________________
    ### density plot at bottom of panel if rug = "density" in function call

    if (rug == "density") {

        # calculate the coordinates of the density function

        density <- density( data[,term.label] )

        # the height of the densityity curve

        max.density <- max(density$y)

        # Get the boundaries of the plot to
        # put the density polygon at the x-line

        plot_coordinates <- par("usr")

        # get the "length" and range of the y-axis
        y.scale <- plot_coordinates[4] - plot_coordinates[3]

        # transform the y-coordinates of the density
        # to the lower 10% of the plotting panel

        density$y <- (0.1 * y.scale / max.density) * density$y + plot_coordinates[3]

        ## plot the polygon

        polygon( density$x , density$y , border = F , col = col.dens) }


    # plot white lines (background color) for 2.5%tile, 1Q, median, 3Q and 97.5%tile through confidence shade and density plot
    #  axis( side = 1 , at = quantiles , labels = FALSE , lwd = 0 , col.ticks = "white"  , lwd.ticks = 2 , tck = 1 )
    # plot white line ( background color for HR = 1 reference
    if (ref) {
        axis( side = 2 , at = 0 , labels = FALSE , lwd = 0 , col.ticks = col.ref  , lwd.ticks = lwd.ref , lty.ticks = lty.ref , tck = 1 )}


    ### _______________ rug = "ticks" ____________________________________________________
    ### rug plot if "ticks" is specified in function call
    if (rug == "ticks") {
        # rugs at datapoints
        axis(side = 1 , line = -1.2 , at = jitter(data[,term.label]) , labels = F , tick = T , tcl = 0.8 , lwd.ticks = 0.1 , lwd = 0)
        # rugs and labels at 1Q, median and 3Q
        axis(side = 1 , line = -1.0 , at = fivenum(data[,term.label])[2:4], lwd = 0 , tick = T, tcl = 1.2 , lwd.ticks = 1 , col.ticks = "black" , labels = c("Quartile 1","Median","Quartile 3"), cex.axis = 0.7, col.axis = "black" , padj = -2.8)
        axis(side = 1 , line = 0.0 , at = fivenum(data[,term.label])[2:4], lwd = 0 , tick = T, tcl = 0.2 , lwd.ticks = 1 , col.ticks = "black", labels = FALSE)
    }


    # last but not least the main plotting line for the smooth estimate:
    # ___________ main plot _____________
    lines(data[,term.label][i.forw] , data$fit[i.forw], col = col.term, lwd = lwd.term)
    # ___________ main plot _____________


    # plot the axes
    if (axes){
        axis(side = 1)
        axis(side = 2 , at = axTicks(2) , label = round( exp(axTicks(2)) , digits = digits))
    }


    # plot a box around plotting panel if specified - not plotted by default
    box(bty = bty)
}


# pacman::p_load(rstpm2) # for brcancer data
# pacman::p_load(TH.data) # same as brcancer data but with proper names
# pacman::p_load(splines)
# devtools::load_all("~/git_repositories/casebase/")
# pacman::p_load(visreg)
#
# # data("brcancer")
# # factoring hormon, breaks the plot method for rstpm2 with interactions!!!
# # brcancer <- transform(brcancer, hormon=factor(hormon, levels = 0:1, labels = c("no","yes")))
# data("GBSG2")
#
# GBSG2 <- transform(GBSG2, hormon=as.numeric(GBSG2$horTh)-1)
# GBSG2 <- transform(GBSG2, meno=as.numeric(GBSG2$menostat)-1)
#
#
# mod_cb_tvc <- fitSmoothHazard(cens ~ hormon*nsx(log(time), df = 3),
#                               data = GBSG2,
#                               time = "time")
#
#
# ci <- function(hazard.params, data) {
#
#     if (is.null(hazard.params[["by"]]))
#         stop("'by' argument needs to be specified in hazard.params argument when type='hr'")
#
#     # name of by argument
#     by_arg <- hazard.params[["by"]]
#
#     # name of xaxis variable
#     x_arg <- hazard.params[["xvar"]]
#
#     if (length(x_arg) > 1) {
#
#         warning("more than one xvar supplied. Only plotting hazard ratio for first element.")
#         x_arg <- x_arg[1]
#
#     }
#
#     ind <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
#     d <- data[ind,]
#
#     fit <- fitSmoothHazard(cens ~ hormon*nsx(log(time), df = 3),
#                            data = d,
#                            time = "time")
#
#     # actual by argument data
#     by_data_vector <- fit[["data"]][[by_arg]]
#
#     # need to check what happens if by is character vector
#     by_unique_values <- sort(unique(by_data_vector))
#
#     tt <- do.call("visreg", utils::modifyList(
#         list(fit = fit,
#              trans = exp,
#              plot = FALSE,
#              rug = FALSE,
#              alpha = 1,
#              partial = FALSE),
#         hazard.params))
#
#     # second level
#     t1 <- tt[["fit"]][which(tt[["fit"]][[by_arg]]==by_unique_values[[2]]),]
#     t1_visreg <- t1[order(t1[[x_arg]]),]
#
#     # reference level
#     t0 <- tt[["fit"]][which(tt[["fit"]][[by_arg]]==by_unique_values[[1]]),]
#     t0_visreg <- t0[order(t1[[x_arg]]),]
#
#
#     hazard_ratio <- cbind(time = t1_visreg[[x_arg]],
#                           HR = t1_visreg[["visregFit"]] / t0_visreg[["visregFit"]])
#
#     return(hazard_ratio)
# }
#
# system.time({
# res <- lapply(1:10, function(i) ci(hazard.params = list(xvar = "time",
#                                                         cond = list(time = quantile(GBSG2$time)),
#                                                         by = "hormon"), data = GBSG2))
# })
#
# pacman::p_load(doParallel)
# doParallel::registerDoParallel(cores = 8)
# system.time({
# res <- foreach(d = 1:100) %dopar%
#     ci(hazard.params = list(xvar = "time",
#                             by = "hormon"), data = GBSG2)
# })
#
# plot(res[[1]])
#
# mod_cb_tvc <- fitSmoothHazard(cens ~ hormon*nsx(log(time), df = 3),
#                               data = GBSG2,
#                               time = "time")
# # mod_rstpm2_tvc <- rstpm2::stpm2(Surv(time,cens == 1) ~ hormon,
# #                                 data = GBSG2, df = 3,
# #                                 tvc=list(hormon = 3))
#
# summary(mod_cb_tvc)
# class(mod_cb_tvc)
#
# mod_cb_tvc$timeVar
#
# dev.off()
# par(mfrow = c(2,1))
# plot(mod_cb_tvc, type = "hazard",
#      hazard.params = list(xvar = c("time"),
#                           by = "hormon",
#                           alpha = 0.05,
#                           ylab = "Hazard"))
#
# plot(mod_cb_tvc, type = "hr",
#      hazard.params = list(xvar = "time",
#                           by = "hormon",
#                           alpha = 0.05,
#                           ylab = "Hazard Ratio"))
#
# DT <- do.call(rbind, res)
# plot(res[[1]][,"time"], res[[2]][,"time"])
# abline(a=0,b=1)
# all.equal(res[[1]][,"time"], res[[2]][,"time"])
# all.equal(res[[1]][,"time"], res[[20]][,"time"])
# lapply(res, lines, col = "grey")
#
#
#
# bs.out <- boot::boot(data=GBSG2,
#                      statistic=ci,
#                      R=3,
#                      hazard.params = list(xvar = "time",
#                                           by = "hormon"))
#
# bs.out
#
#
# class(mod_cb_tvc)
#
# library(car)
#
# set.seed(1234)
# duncan.boot <- car:::Boot.default(mod_cb_tvc, R = 19)
#
# View(car:::Boot.default)
# mod_cb_tvc$data
#
# tt <- summary(duncan.boot, high.moments=TRUE)
#
# library(magrittr)
# library(knitr)
#
#
# hist(duncan.boot, legend="separate")
#
#
#
# termss <- predict(mod_cb_tvc , type="terms" , se.fit = TRUE, newdata = data.frame(hormon = 0, time = quantile(mod_cb_tvc$data$time), offset = 0))
# termss$fit %>% head
#
# par(mfrow = c(2,1))
# plot(cbind(quantile(mod_cb_tvc$data$time), exp(termss$fit[,1])), type = "l")
#
# exp(rowSums(termss$fit) + attr(termss$fit,"constant"))
#
# tpp <- plot(mod_cb_tvc, type = "hr",
#      hazard.params = list(xvar = c("time"),
#                           by = "hormon",
#                           alpha = 0.05,
#                           ylab = "Hazard"))
#
# tpp
# summary(mod_cb_tvc)
# exp(-coef(mod_cb_tvc))
#
#
# library(survival)
# survobj<-Surv(time, cens)
#
# str(GBSG2)
#
# model <- coxph(Surv(time, cens)~age+hormon, data=GBSG2)
# summary(model)
#
# model <- glm(cens ~ nsx(log(time),df=3)*hormon, data = GBSG2)
# termss <- predict(model , type="terms" , se.fit = TRUE)
# termss$fit %>% head
#
# plotHR(model, terms="time", se = TRUE , rug = "density" ,
#        xlab = "" , ylab = "" ,
#        # xlim=c(), ylim=c(),
#        col.term = "#08519C",
#        lwd.term = "2",
#        col.se = "#DEEBF7",
#        cex = "1.5" ,
#        bty = "L" ,
#        axes = TRUE)
#
#
# test.data = data.frame(y = c(0,0,0,1,1,1,1,1,1), x=c(1,2,3,1,2,2,3,3,3))
# model = glm(y~(x==1)+(x==2), family = 'binomial', data = test.data)
# summary(model)
# X <- model.matrix(y~(x==1)+(x==2), data = test.data)
#
# beta <- model$coefficients
# beta * t(X)
#
# model$coefficients
# tt <- predict(model, type = "terms", se.fit = TRUE)
# rowSums(tt$fit) + attr(tt$fit, "constant")
# predict(model, type = "link", se.fit = TRUE)
#
#
# plot(ppp$coefficients[,"Std. Error"],
#      tt$bootSE)
# abline(a=0,b=1)

