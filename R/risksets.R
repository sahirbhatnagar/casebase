
# Simulate censored survival data for two outcome types from Weibull distributions:
library(data.table)
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
DT <- data.table(z=rbinom(nobs, 1, 0.5))
DT[,`:=` ("t_event"=rweibull(nobs, a1, b1 * exp(z * c1)^(-1/a1)),
          "t_comp"=rweibull(nobs, a2, b2 * exp(z * c2)^(-1/a2)),
          "tlim"=10)]
DT[,`:=`("e"=1 * (t_event < t_comp) + 2 * (t_event >= t_comp),
         "t"=pmin(t_event, t_comp))]
DT[t >= tlim, e:=0]
DT[t >= tlim, t:=tlim]
DT[, idx:=order(t)]
# put the e=0 people at the bottom of the population-time plot and people with
# short values of t at the top
DT[ DT[["idx"]], ycoord:=(nobs:1)]
# sample y coordinates for each event, so that we can see the incidence density
# on population-time plots. Sampling from people with e=0 or 2 who have an
# observed time t greater than that of a fixed individual who had the event
DT[e==1, yc:=lapply(t, function(i){sample(DT[t>=i & e!=1,ycoord],1)})]


# get row number of subjects who have an event==1, and either covariate value
cond <- DT[e==1 & z %in% c(0,1), which=T]

# get event time of subjects who have an event==1, and either covariate value
etimes <- DT[e==1 & z %in% c(0,1), t]


# DT[e==1, ycoordnew:=mapply(function(x,y) {DT[ycoord >= yc & ycoord < i,ycoord ]+1}  )]
# mapply
#     ycoord[DT[["t"]]>=i & DT[["e"]]!=1],1)}  )]
#
# for (i in 1:nobs) {
#     if (e[i] == 1) {
#         # yc: for a given subject i with event=1, sample the id numbers from
#         # subjects who have an event time (either from being censored or event
#         # type=2) greater than subject i
#         yc <- sample(ycoord[t >= t[i] & e != 1], 1)
#
#         ycoord[ycoord >= yc & ycoord < ycoord[i]] <- ycoord[ycoord >= yc & ycoord < ycoord[i]] + 1
#         ycoord[i] <- yc
#
#     }
# }

DT[,table(e)]

op <- par(mar=c(4.5,4.5,1,1), lwd=2)
with(DT,plot(1,1,type='n', xlim=c(0,tlim[1]), ylim=c(1,nobs), xlab='Follow-up years', ylab='Population'))
with(DT,segments(rep(0.0, nobs), ycoord, t, ycoord, col='gray'))
#with(DT[z==1],segments(rep(0.0, nobs), ycoord, t, ycoord, col='yellow'))
#with(DT[z==0],points(t[e==1], yc[e==1], pch=20, col='blue', cex=0.5))
with(DT,points(t[e==1], yc[e==1], pch=20, col='blue', cex=0.5))
par(op)


set.seed(1)
nobs <- 5000
# binary covariate
z <- rbinom(nobs, 1, 0.5)

# simulated times based on binary covariate
t1 <- rweibull(nobs, a1, b1 * exp(z * c1)^(-1/a1))
t2 <- rweibull(nobs, a2, b2 * exp(z * c2)^(-1/a2))

# censoring time
tlim <- 10

# event rate
mean(t1 < tlim)
mean(t2 < tlim)

# event type 0-censored, 1-event of interest, 2-competing event
e <- 1 * (t1 < t2) + 2 * (t1 >= t2)

# which event happens first
t <- pmin(t1, t2)

# censor those who exceed study period
e[t >= tlim] <- 0
t[t >= tlim] <- tlim
table(e)


# Reorder for plotting:
idx <- order(t)
ycoord <- rep(NA, nobs)
for (i in 1:nobs) {
    ycoord[idx[i]] <- (nobs:1)[i]
}

idx[1:10]
ycoord[1:10]


for (i in 1:nobs) {
    if (e[i] == 1) {
        # yc: for a given subject i with event=1, sample the id numbers from
        # subjects who have an event time (either from being censored or event
        # type=2) greater than subject i
        yc <- sample(ycoord[t >= t[i] & e != 1], 1)
        ycoord[ycoord >= yc & ycoord < ycoord[i]] <- ycoord[ycoord >= yc & ycoord < ycoord[i]] + 1
        ycoord[i] <- yc

    }
}


# Population-time plots:

pdf(file.path(getwd(), 'riskset0.pdf'), width=7, height=6, paper='special')
op <- par(mar=c(4.5,4.5,1,1), lwd=2)
plot(1,1,type='n', xlim=c(0,tlim), ylim=c(1,nobs), xlab='Follow-up years', ylab='Population')
segments(0.0, 1, 0.0, nobs, col='gray')
par(op)
dev.off()

pdf(file.path(getwd(), 'riskset1.pdf'), width=7, height=6, paper='special')
op <- par(mar=c(4.5,4.5,1,1), lwd=2)
plot(1,1,type='n', xlim=c(0,tlim), ylim=c(1,nobs), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, nobs), ycoord, t, ycoord, col='gray')
par(op)
dev.off()

pdf(file.path(getwd(), 'riskset2.pdf'), width=7, height=6, paper='special')
op <- par(mar=c(4.5,4.5,1,1), lwd=2)
plot(1,1,type='n', xlim=c(0,tlim), ylim=c(1,nobs), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, nobs), ycoord, t, ycoord, col='gray')
points(t[e==1], ycoord[e==1], pch=20, col='blue', cex=0.5)
par(op)
dev.off()

last <- 80
idx <- ycoord <= last
pdf(file.path(getwd(), 'riskset3.pdf'), width=7, height=6, paper='special')
op <- par(mar=c(4.5,4.5,1,1), lwd=2)
plot(1,1,type='n', xlim=c(0,tlim), ylim=c(1,last), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, nobs)[idx], ycoord[idx], t[idx], ycoord[idx], col='gray')
points(t[e==1], ycoord[e==1], pch=20, col='blue', cex=0.5)
par(op)
dev.off()

pdf(file.path(getwd(), 'riskset4.pdf'), width=7, height=6, paper='special')
op <- par(mar=c(4.5,4.5,1,1), lwd=2)
plot(1,1,type='n', xlim=c(0,tlim), ylim=c(1,last), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, nobs)[idx], ycoord[idx], t[idx], ycoord[idx], col='gray')
points(t[e==1], ycoord[e==1], pch=20, col='blue', cex=0.5)
for (ti in (t[idx])[e[idx]==1]) {
    points(rep(ti, last)[t[idx] > ti], (ycoord[idx])[t[idx] > ti], pch=20, col='gray', cex=0.5)
}
par(op)
dev.off()

pdf(file.path(getwd(), 'riskset5.pdf'), width=7, height=6, paper='special')
op <- par(mar=c(4.5,4.5,1,1), lwd=2)
plot(1,1,type='n', xlim=c(0,tlim), ylim=c(1,last), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, nobs)[idx], ycoord[idx], t[idx], ycoord[idx], col='gray')
points(t[e==1], ycoord[e==1], pch=20, col='blue', cex=0.5)
for (ti in (t[idx])[e[idx]==1]) {
    points(rep(ti, last)[t[idx] > ti], (ycoord[idx])[t[idx] > ti], pch=20, col='gray', cex=0.5)
}
for (ti in (t[idx])[e[idx]==1]) {
    sampled <- sample((1:last)[t[idx] > ti], 1)
    points(rep(ti, last)[sampled], (ycoord[idx])[sampled], pch=20, col='red', cex=0.5)
}
par(op)
dev.off()

pdf(file.path(getwd(), 'riskset6.pdf'), width=7, height=6, paper='special')
op <- par(mar=c(4.5,4.5,1,1), lwd=2)
plot(1,1,type='n', xlim=c(0,tlim), ylim=c(1,last), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, nobs)[idx], ycoord[idx], t[idx], ycoord[idx], col='gray')
points(t[e==1], ycoord[e==1], pch=20, col='blue', cex=0.5)
for (ti in (t[idx])[e[idx]==1]) {
    points(rep(ti, last)[t[idx] > ti], (ycoord[idx])[t[idx] > ti], pch=20, col='gray', cex=0.5)
}
for (ti in (t[idx])[e[idx]==1]) {
    sampled <- sample((1:last)[t[idx] > ti], 5)
    points(rep(ti, last)[sampled], (ycoord[idx])[sampled], pch=20, col='red', cex=0.5)
}
par(op)
dev.off()

pdf(file.path(getwd(), 'riskset7.pdf'), width=7, height=6, paper='special')
op <- par(mar=c(4.5,4.5,1,1), lwd=2)
plot(1,1,type='n', xlim=c(0,tlim), ylim=c(1,last), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, nobs)[idx], ycoord[idx], t[idx], ycoord[idx], col='gray')
points(t[e==1], ycoord[e==1], pch=20, col='blue', cex=0.5)
for (ti in (t[idx])[e[idx]==1]) {
    points(rep(ti, last)[t[idx] > ti], (ycoord[idx])[t[idx] > ti], pch=20, col='red', cex=0.5)
}
par(op)
dev.off()



