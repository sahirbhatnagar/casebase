library(data.table)
library(magrittr)
library(ggplot2)

rm(list=ls())

popTimeData <- popTime(data = veteran)


bmt <- read.csv("https://raw.githubusercontent.com/sahirbhatnagar/casebase/master/inst/extdata/bmtcrr.csv")

head(bmt)
str(bmt)
popTimeData <- popTime(data = bmt, time = "ftime", event = "Phase")
str(popTimeData)
popTimeData[, table(`event status`)]
popTimeData[, table(event)]

popTimeData <- popTime(data = DT, time = "time")



data("heart")
str(heart)
heart$time <- heart$stop - heart$start
str(heart)
popTimeData <- popTime(data = heart, exposure = "transplant")

data("stanford2")
str(stanford2)

data("cancer")
str(cancer)




cancer$status <- cancer$status - 1
# error because not enough cases in subcategories
popTimeData <- popTime(data = cancer, exposure = "ph.ecog")

popTimeData <- popTime(data = cancer, exposure = "sex")
xtabs(~ph.ecog+status, cancer)


popTimeData <- popTime(data = DTsim, time = "time", event = "event", exposure = "z")
popTimeData <- popTime(data = bmt, time = "ftime", exposure = "Sex")




popTimeData <- popTime(data = stanford2)
class(popTimeData)
plot(popTimeData)

popTimeData



p1 <- ggplot(popTimeData, aes(x=0, xend=time, y=ycoord, yend=ycoord)) +
    geom_segment(size=.5, colour="grey80") +
    xlab("Follow-up years") +
    ylab("Population") +
    theme_bw() +
    theme(axis.text=element_text(size=12, face='bold'), legend.position = "bottom",
          legend.title=element_blank())

p1 + facet_grid(~transplant)
p1 + geom_point(aes(x=time, y=yc, colour = `event status`), data = popTimeData[event==1],
                size=1) +
    facet_grid(~transplant)
#facet_wrap(~transplant, ncol = 1)

p1 + geom_point(aes(x=time, y=yc), data = popTimeData[event==1],
                size=1, color = "red") +
    facet_grid(~transplant)



# Simulated Data Example --------------------------------------------------

set.seed(1)
nobs <- 5000

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
DTsim <- data.table(ID = seq_len(nobs), z=rbinom(nobs, 1, 0.5))

setkey(DTsim, ID)

DTsim[,`:=` (event_time = rweibull(nobs, a1, b1 * exp(z * c1)^(-1/a1)),
             competing_time = rweibull(nobs, a2, b2 * exp(z * c2)^(-1/a2)),
             end_of_study_time = eost)]

DTsim[,`:=`(event = 1 * (event_time < competing_time) + 2 * (event_time >= competing_time),
            time = pmin(event_time, competing_time))]

DTsim[time >= end_of_study_time, event := 0]

DTsim[time >= end_of_study_time, time:=end_of_study_time]
