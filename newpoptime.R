devtools::load_all()
library(splines)
library(visreg)

# support data ------------------------------------------------------------
data("support")

mod_cb <- casebase::fitSmoothHazard(death ~ bs(d.time)*sex,
                                    time = "d.time",
                                    event = "death",
                                    data = support,
                                    ratio = 10)
summary(mod_cb)


hazardPlot(mod_cb, newdata = data.frame(sex = "male"), ci = F)
hazardPlot(mod_cb, newdata = data.frame(sex = "female"), add = T)

absoluteRisk()

visreg(mod_cb,
       xvar = "d.time",
       by = "sex",
       scale = "response",
       overlay = T,
       # data = mod_cb$data,
       plot = T)

visreg(mod_cb,
       xvar = "sex",
       by = "d.time",
       type = "contrast",
       scale = "response",
       data = mod_cb$data,
       plot = T)

mod_cb$data$offset <- 0
visreg(mod_cb,
       xvar = "d.time",
       by = "sex",
       # scale = "linear",
       scale = "response",
       data = mod_cb$data,
       plot = T)

tt$fit
tt$res
# xtrans = exp)






# simdat ------------------------------------------------------------------

devtools::load_all()

library(casebase)
library(splines)
library(visreg)
data("simdat")
mod_cb <- casebase::fitSmoothHazard(status ~ trt + ns(I(log(eventtime)), df = 3) +
                                        trt:ns(I(log(eventtime)),df=1),
                                    time = "eventtime",
                                    data = simdat,
                                    ratio = 10,
                                    family = "glm")

results0 <- hazardPlot(object = mod_cb, newdata = data.frame(trt = 0),
                       ci.lvl = 0.95, ci = TRUE, lty = 1, line.col = 1, lwd = 2)
head(results0)
hazardPlot(object = mod_cb, newdata = data.frame(trt = 1), ci = TRUE,
           ci.lvl = 0.95, add = TRUE, lty = 2, line.col = 2, lwd = 2)
legend("topleft", c("trt=0","trt=1"),lty=1:2,col=1:2,bty="y", lwd = 2)

mod_cb$data$offset <- 0
visreg(mod_cb,
       xvar = "eventtime",
       by = "trt",
       scale = "response",
       data = mod_cb$data,
       plot = T)


data("bmtcrr")
head(bmtcrr)
str(bmtcrr)

data("support")
?support
popTimeData <- popTime(data = bmtcrr, time = "ftime", event = "Status", exposure = "D")
# popTimeData <- popTime(data = bmtcrr, time = "ftime", event = "Status")
data("ERSPC")
head(ERSPC)

library(splines)
devtools::load_all()
data("support")
mod_cb <- casebase::fitSmoothHazard(death ~ bs(d.time)+sex,
                                    time = "d.time",
                                    event = "death",
                                    data = support,
                                    ratio = 10)
summary(mod_cb)
support$sex
hazardPlot(mod_cb, newdata = data.frame(sex = "male"), ci = F)
hazardPlot(mod_cb, newdata = data.frame(sex = "female"), add = T)


pacman::p_load(effects)
ggeffects::ggpredict(mod_cb)
str(mod_cb)
head(mod_cb$data)
head(mod_cb$originalData)




# load packages and data ---------------------------------------------------

devtools::load_all("~/git_repositories/casebase/")
pacman::p_load(splines)
pacman::p_load(ggeffects)
pacman::p_load(visreg)
data("simdat")
head(simdat)

# Fit the model on simdat from casebase -----------------------------------

mod_cb <- casebase::fitSmoothHazard(status ~ bs(eventtime,degree = 8) * trt,
                                    time = "eventtime",
                                    event = "status",
                                    data = simdat,
                                    ratio = 10)
summary(mod_cb)
plot(mod_cb)
str(mod_cb)

mod_cb

# using casebase hazardPlot ------------------------------------------------

t1_hazardplot <- hazardPlot(object = mod_cb, newdata = data.frame(trt = 1),
                       ci.lvl = 0.95, ci = TRUE, lty = 1, line.col = 1, lwd = 2)

t0_hazardplot <- hazardPlot(object = mod_cb, newdata = data.frame(trt = 0), ci = TRUE,
           ci.lvl = 0.95, add = TRUE, lty = 2, line.col = 2, lwd = 2)
legend("topleft", c("trt=1","trt=0"),lty=1:2,col=1:2,bty="y", lwd = 2)


# using visreg ------------------------------------------------------------

mod_cb$data$offset <- 0
tt <- visreg(mod_cb,
       xvar = "eventtime",
       by = "trt",
       scale = "linear",
       trans = exp,
       data = mod_cb$data,
       plot = T,
       overlay = TRUE)

# get the predicted hazards by trt status and order the
# times so you can plot them
t1 <- tt$fit[which(tt$fit$trt==1),]
t0 <- tt$fit[which(tt$fit$trt==0),]
t1_visreg <- t1[order(t1$eventtime),]
t0_visreg <- t0[order(t1$eventtime),]



# ggpredict ---------------------------------------------

# this stops predictions at eventtimes <0.5
dat1 <- ggpredict(mod_cb, terms = c("eventtime [all]","trt","offset [0]"))
dat1$x %>% range

# this works, but produces slightly different results than visreg or hazardPlot
dat <- ggpredict(mod_cb, terms = c("eventtime [0:5 by=.1]","trt","offset [0]"))
plot(dat)

# get the predicted hazards by trt status. they are already in order
# wih ggpredict
t1_ggpredict <- dat$predicted[which(dat$group==1)]
t0_ggpredict <- dat$predicted[which(dat$group==0)]



# compare hazard ratios for treatment -----------------------------------------

plot(t0_hazardplot$eventtime, t1_hazardplot$predictedhazard / t0_hazardplot$predictedhazard,
     type = "l", col = 1, lty = 1, ylab = "hazard ratio for trt=1 vs. trt=0", lwd = 2, ylim = c(0,1.1))
abline(a=1, b=0, col = "grey80")
lines(t0_visreg$eventtime, t1_visreg$visregFit / t0_visreg$visregFit,
      type = "l", col = 2, lty = 2, lwd = 2)
lines(unique(dat$x), t1_ggpredict / t0_ggpredict,
      type = "l", col = 3, lty = 3, lwd = 2)

legend("bottomright", c("hazardPlot","visreg","ggpredict"),
       lty=1:3,col=1:3,bty="y", lwd = 2)







summary(mod_cb)

t1$eventtime == t0$eventtime

plot()

mod_cb$data <- mod_cb$originalData

plot(dat)

mod_cb$data %>% colnames()
mod_cb$originalData %>% colnames()
mod_cb$data <- mod_cb$originalData



pacman::p_load(visreg)

mod_cb$originalData$offset = mod_cb$offset[1]
visreg(mod_cb,
       xvar = "d.time",
       by = "sex",
       scale = "linear",
       data = mod_cb$originalData)#,
       # xtrans = exp)


dev.off()
modcb$originalData %>% str
modcb$data <- modcb$originalData
modcb$data %>% str
modcb$originalData %>% str
plot(dat, facet = TRUE)

head(modcb$data)
head(modcb$originalData)
library(survival)
data(veteran)
head(veteran)
veteran$status %>% table

modcb <- fitSmoothHazard(status ~ ., data = veteran)
summary(modcb)

data("simdat")
library(splines)
mod_cb <- casebase::fitSmoothHazard(status ~ trt + ns(log(eventtime), df = 3) +
                                        trt:ns(log(eventtime),df=1),
                                    time = "eventtime",
                                    data = simdat,
                                    ratio = 100,
                                    family = "glm")

results0 <- hazardPlot(object = mod_cb, newdata = data.frame(trt = 0),
                       ci.lvl = 0.95, ci = TRUE, lty = 1, line.col = 1, lwd = 2)
head(results0)
hazardPlot(object = mod_cb, newdata = data.frame(trt = 1), ci = TRUE,
           ci.lvl = 0.95, add = TRUE, lty = 2, line.col = 2, lwd = 2)
legend("topleft", c("trt=0","trt=1"),lty=1:2,col=1:2,bty="y", lwd = 2)


summary(modcb)
hazardPlot(modcb, newdata = ERSPC)
?fitSmoothHazard
?hazardPlot
attr(popTimeData, "exposure")
attr(popTimeData, "call")

plot(popTimeData,
     add.base.series = T,
     comprisk = T,
     add.competing.event = T,
     ratio = 1,
     legend = TRUE,
     facet.params = list(ncol=1),
     theme.params = list(legend.position = "top"),
     casebase.theme = F
)

str(popTimeData)






library(casebase)
library(ggplot2)

ERSPC$DeadOfPrCa <- factor(ERSPC$DeadOfPrCa,
                           levels = c(0, 1),
                           labels = c("Censored", "PrCa Death"))
popTimeData <- popTime(data = ERSPC,
                       time = "Follow.Up.Time",
                       event = "DeadOfPrCa")
head(popTimeData)
p1 <- ggplot()

# doesnt work
nogo <- list(data = popTimeData[event == 1], mapping = aes(x = time, y = yc, colour = `event status`), colour = "green")
p1 + do.call("geom_point", nogo)

# works
works <- list(data = popTimeData[event == 1], mapping = aes(x = time, y = yc, colour = `event status`))
p1 + do.call("geom_point", works) + scale_color_manual(values = c("PrCa Death" = "blue"))


dim(popTimeData)
sampleCaseBase(popTimeData)
sampleData <- ERSPC %>%
    arrange(desc(Follow.Up.Time)) %>%
    mutate(id = row_number()) %>%
    sampleCaseBase(time = "Follow.Up.Time",
                   event = "DeadOfPrCa",
                   comprisk = FALSE,
                   ratio = 1)
plot(popTimeData) + geom_point(aes(x = Follow.Up.Time, y = id),
                               data = filter(sampleData,
                                             DeadOfPrCa == 0),
                               colour = "black",
                               size = 0.5, inherit.aes = FALSE)


huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
h <- ggplot(huron, aes(year))

h + geom_ribbon(aes(ymin=0, ymax=level))


popTimeData %>% head

huron %>% head


library(casebase)
library(ggplot2)


devtools::load_all()
ERSPC$ScrArm <- factor(ERSPC$ScrArm,
                       levels = c(0, 1),
                       labels = c("Control group", "Screening group"))
ERSPC$DeadOfPrCa <- factor(ERSPC$DeadOfPrCa,
                       levels = c(0, 1),
                       labels = c("Censored", "PrCa Death"))
popTimeData <- popTime(data = ERSPC,
                       time = "Follow.Up.Time",
                       event = "DeadOfPrCa")

plot(popTimeData)

devtools::load_all()
cols <- c("Case series" = "#D55E00", "Competing event" = "#009E73", "Base series" = "#0072B2")
plot(popTimeData,
     casebase.theme = TRUE,
     add.case.series = TRUE,
     ratio = 1,
     add.base.series = TRUE,
     legend = TRUE,
     base.params = list(show.legend = TRUE))#,
     legend.params = list(name = element_blank(),
                     breaks = c("Case series", "Competing event", "Base series"),
                     values = cols))
     # case.params = list(mapping = aes(x = time, y = yc, colour = "Relapse")),
     # base.params = list(mapping = aes(x = time, y = ycoord, colour = "controls")),
     # legend.params = list(breaks = c("Relapse", "controls"), values = c("Relapse" = "green","controls" = "red")))


popTimeData <- popTime(data = bmtcrr, time = "ftime", event = "Status")
cols <- c("Case series" = "black", "Competing event" = "#009E73", "Base series" = "#0072B2")
p1 <- plot(popTimeData,
     casebase.theme = TRUE,
     add.case.series = TRUE,
     ratio = 1,
     add.base.series = TRUE,
     legend = TRUE,
     comprisk = TRUE,
     add.competing.event = FALSE,
     legend.params = list(name = element_blank(),
                     breaks = c("Case series", "Competing event", "Base series"),
                     values = cols),
     theme.params = list(legend.position = "none"))

dev.off()
p1

p1 + labs(caption = c("hellow worls"), title = "Poptime", x = "nine")
     # case.params = list(mapping = aes(x = time, y = yc, colour = "Relapse")),
     # base.params = list(mapping = aes(x = time, y = ycoord, colour = "controls")),
     # legend.params = list(breaks = c("Relapse", "controls"), values = c("Relapse" = "green","controls" = "red")))


h <- ggplot(popTimeData)


devtools::load_all()
# library(casebase)
popTimeData <- popTime(data = bmtcrr, time = "ftime", event = "Status")
plot(popTimeData,
     add.case.series = TRUE,
     add.base.series = TRUE,
     add.competing.event = FALSE,
     comprisk = TRUE,
     # case.params = list(show.legend = TRUE),
     ratio = 1,
     legend = T,
     casebase.theme = TRUE,
     theme.params = list(legend.position = "bottom"))



bt <- data.table::as.data.table(bmtcrr)
popTimeData <- popTime(data = bt, time = "ftime", event = "Status")



head(popTimeData)


rlang::last_error()
?popTime

head(popTimeData)
popTimeData$comprisk.event <- 0
popTimeData$comprisk.event[popTimeData$event == 1] <- 2
popTimeData$comprisk.event[popTimeData$event == 2] <- 1

table(popTimeData$comprisk.event, popTimeData$event)

str(popTimeData)
popTime(data = popTimeData, time = "time", event = "comprisk.event")




# Ribbon ------------------------------------------------------------------

h + geom_ribbon(aes(x = time, ymin=0, ymax=ycoord)) +
    theme_minimal()

pdf("geom_ribbon.pdf")
h + geom_ribbon(aes(x = time, ymin=0, ymax=ycoord)) +
    theme_minimal()
dev.off()



# Segment -----------------------------------------------------------------

h + geom_segment(aes(x = 0, xend = time, y = ycoord, yend = ycoord)) +
    theme_minimal()

pdf("geom_segment.pdf")
h + geom_segment(aes(x = 0, xend = time, y = ycoord, yend = ycoord)) +
    theme_minimal()
dev.off()



# File size in kB ---------------------------------------------------------

file.info("geom_ribbon.pdf")$size/1000
file.info("geom_segment.pdf")$size/1000



x_var <- "displ"
aes(x_var)


aes_(quote(displ), quote(ty), quote(uu))
#> Aesthetic mapping:
#> * `x` -> `displ`
aes_(as.name(x_var))
#> Aesthetic mapping:
#> * `x` -> `displ`
aes_(parse(text = x_var)[[1]])
#> Aesthetic mapping:
#> * `x` -> `displ`

f <- function(x_var) {
    aes_(substitute(x_var))
}
f(displ)
#> Aesthetic mapping:
#> * `x` -> `displ`



library(ggplot2)

mtcars$wt2 <- mtcars$wt*0.9
mtcars$mpg2 <- mtcars$mpg*0.9

mp <- list(size = 3, color = "red", pch = 10, aes(x = wt2, y = mpg2))

add <- FALSE

tt <- list(
    if (add)
        do.call("geom_point",mp)
)

p <- ggplot(mtcars, aes(wt, mpg))
p + geom_point() + tt



library(magrittr)
do.call("geom_point",mp)


cols <- c("LINE1"="#f04546","LINE2"="#3591d1","BAR"="#62c76b")
a <-c("S1","S2","S3","S4","S5","S6","S7","S8","S9") #names
b <-c(0.23,0.26,0.55,0.56,0.36,0.23,0.18,0.06,0.04) #mean t0
c <-c(0.64,0.6,0.81,1.4,0.89,0.55,0.48,0.22,0.09) #mean t1
d <-c(0.20,0.23,0.52,0.53,0.33,0.20,0.15,0.04,0.03) #SD low t0
e <-c(0.26,0.29,0.58,.59,0.39,0.26,0.21,0.08,0.05) #SD high t0
f <-c(0.67,0.63,0.86,1.44,0.93,0.59,0.51,0.25,0.10) #SD high t1
g <-c(0.61,0.57,0.78,1.36,0.85,0.53,0.45,0.19,0.08) #SD low t1
h <-c(0.41,0.34,0.26,0.84,0.53,0.32,0.30,0.16,0.05) #absolute change

data <- data.frame(a,b,c,d,e,f,g,h)
ggplot(data=data,aes(x=a)) +
    geom_bar(stat="identity", aes(y=h, fill = "BAR"),colour="#333333")+ #green
    geom_line(aes(y=b,group=1, colour="LINE1"),size=1.0) +   #red
    geom_point(aes(y=b, colour="LINE1"),size=3) +           #red
    geom_errorbar(aes(ymin=d, ymax=e, colour="LINE1"), width=0.1, size=.8) +
    geom_line(aes(y=c,group=1,colour="LINE2"),size=1.0) +   #blue
    geom_point(aes(y=c,colour="LINE2"),size=3) +           #blue
    geom_errorbar(aes(ymin=f, ymax=g,colour="LINE2"), width=0.1, size=.8) +
    scale_colour_manual(name="Error Bars",values=cols) + scale_fill_manual(name="Bar",values=cols) +
    ylab("Symptom severity") + xlab("PHQ-9 symptoms") +
    ylim(0,1.6) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust=-.2)) +
    theme(axis.title.y = element_text(size = 15, vjust=0.3))

set.seed(7)
df <- data.frame(
    x = 1:20,
    y.bar = rpois(20, lambda = 5),
    y.line = rpois(20, lambda = 10)
)
ggplot(data = df,
       mapping = aes(x = x)) +

    # specify fill for bar / color for line inside aes(); you can use
    # whatever label you wish to appear in the legend
    geom_col(aes(y = y.bar, fill = "bar.label")) +
    geom_line(aes(y = y.line, color = "line.label")) +

    xlab("Month of year") +
    scale_y_continuous(name = "Daily classifications per Spotter") +

    # the labels must match what you specified above
    scale_fill_manual(name = "", values = c("bar.label" = "grey")) +
    scale_color_manual(name = "", values = c("line.label" = "black")) +

    theme_bw()
