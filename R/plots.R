library(data.table)
library(magrittr)

rm(list=ls())
set.seed(1)
nobs <- 5000

# simulation parameters
a1 <- 1.0
b1 <- 200
a2 <- 1.0
b2 <- 50
c1 <- 0.0
c2 <- 0.0

# e event type 0-censored, 1-event of interest, 2-competing event
# t observed time/endpoint
# z is a binary covariate
DT <- data.table(ID = seq_len(nobs), z=rbinom(nobs, 1, 0.5))

setkey(DT, ID)

DT[,`:=` (event_time = rweibull(nobs, a1, b1 * exp(z * c1)^(-1/a1)),
          competing_time = rweibull(nobs, a2, b2 * exp(z * c2)^(-1/a2)),
          end_of_study_time = 10)]
DT[,`:=`(event=1 * (event_time < competing_time) + 2 * (event_time >= competing_time),
         time=pmin(event_time, competing_time))]
DT[time >= end_of_study_time, event := 0]
DT[time >= end_of_study_time, time:=end_of_study_time]


# put the e=0 people at the bottom of the population-time plot and people with
# short values of t at the top
DT[ DT[,order(time)], ycoord:=(nobs:1)]
# sample y coordinates for each event, so that we can see the incidence density
# on population-time plots. Sampling from people with e=0 or 2 who have an
# observed time t greater than that of a fixed individual who had the event
DT[event == 1, yc := lapply(time, function(i) sample(DT[time >= i & event != 1,ycoord],1) )]
DT[ , yc := if(event == 1) lapply(time, function(i) sample(DT[time >= i & event != 1,ycoord],1) ) else NA, by = ID]

for (t in DT)



plot(DT$yc)

DT[, yc:=unlist(yc)]


str(DT)


# get row number of subjects who have an event==1, and either covariate value
cond <- DT[e==1 & z %in% c(0,1), which=T]

# get event time of subjects who have an event==1, and either covariate value
etimes <- DT[e==1 & z %in% c(0,1), t]

