library(data.table)
library(magrittr)
library(ggplot2)
library(survival)

rm(list=ls())

# Veteran Data ------------------------------------------------------------

source("R/plots.R")
source("R/methods.R")

# veteran data in library(survival)
data("veteran")
str(veteran)

# create 'popTime' object
popTimeData <- popTime(data = veteran)

# object of class 'popTime'
class(popTimeData)

# plot method for objects of class 'popTime'
plot(popTimeData)

# stratified by treatment population time plot
veteran <- transform(veteran, trt = factor(trt, levels = 1:2,
                                           labels = c("standard", "test")))

# create 'popTimeExposure' object
popTimeData <- popTime(data = veteran, exposure = "trt")

# object of class 'popTimeExposure'
class(popTimeData)

# plot method for objects of class 'popTimeExposure'
plot(popTimeData)



# Stem Cell Data ----------------------------------------------------------

bmt <- read.csv("https://raw.githubusercontent.com/sahirbhatnagar/casebase/master/inst/extdata/bmtcrr.csv")
str(bmt)

# create 'popTime' object
popTimeData <- popTime(data = bmt, time = "ftime")

# object of class 'popTime'
class(popTimeData)

# plot method for objects of class 'popTime'
plot(popTimeData)

# stratified by Disease population time plot
# Disease (lymphoblastic or myeloblastic leukemia,
# abbreviated as ALL and AML, respectively)

# create 'popTimeExposure' object
popTimeData <- popTime(data = bmt, time = "ftime", exposure = "D")

# object of class 'popTimeExposure'
class(popTimeData)

# plot method for objects of class 'popTimeExposure'
plot(popTimeData)

# stratify by gender
popTimeData <- popTime(data = bmt, time = "ftime", exposure = "Sex")
plot(popTimeData)




# Stanford Heart Transplant Data ------------------------------------------

# data from library(survival)
data("heart")
str(heart)

# create time variable for time in study
heart <- transform(heart,
                   time = stop - start,
                   transplant = factor(transplant,
                                       labels = c("no transplant", "transplant")))

# stratify by transplant indicator
popTimeData <- popTime(data = heart, exposure = "transplant")

# can specify a legend
plot(popTimeData, legend = TRUE)



# NCCTG Lung Cancer Data --------------------------------------------------

# data from library(survival)
data("cancer")
str(cancer)

# since the event indicator 'status' is numeric, it must have
# 0 for censored and 1 for event
cancer <- transform(cancer,
                    status = status - 1,
                    sex = factor(sex, levels = 1:2,
                                 labels = c("Male", "Female")))


# population time plot
# redistributing the red points among those who never experienced an event
# because there are enough available at each time point
popTimeData <- popTime(data = cancer)
plot(popTimeData)

# stratified by sex
popTimeData <- popTime(data = cancer, exposure = "sex")

# can change the plot aesthetics
plot(popTimeData,
     line.width = 0.2, line.colour = "black",
     point.size = 1, point.colour = "cyan")



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
DTsim[,`:=`(event = 1 * (event_time < competing_time) +
                2 * (event_time >= competing_time),
            time = pmin(event_time, competing_time))]
DTsim[time >= end_of_study_time, event := 0]
DTsim[time >= end_of_study_time, time:=end_of_study_time]

# create 'popTime' object
popTimeData <- popTime(data = DTsim, time = "time", event = "event")
plot(popTimeData)

# stratified by binary covariate z
popTimeData <- popTime(data = DTsim, time = "time", event = "event", exposure = "z")

# we can line up the plots side-by-side instead of one on top of the other
plot(popTimeData, ncol = 2)







