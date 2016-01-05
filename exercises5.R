
# 1-3:

library(survival)

data(veteran)
str(veteran)

table(veteran$status)
evtimes <- veteran$time[veteran$status == 1]
hist(evtimes, nclass=30, main='', xlab='Survival time (days)', col='gray90')

hist(evtimes, nclass=30, main='', xlab='Survival time (days)', col='gray90', probability=TRUE)
tgrid <- seq(0, 1000, by=10)
lines(tgrid, dexp(tgrid, rate=1.0/mean(evtimes)), lwd=2, lty='dashed', col='red')

y <- Surv(veteran$time, veteran$status)
karno <- veteran$karno
diagtime <- veteran$diagtime
age <- veteran$age
prior <- as.factor(veteran$prior)
celltype <- factor(as.character(veteran$celltype), levels=c('large', 'squamous', 'smallcell', 'adeno'))
trt <- as.factor(veteran$trt)

library(eha)
model1 <- weibreg(y ~ karno + diagtime + age + prior + celltype + trt, shape=1)
summary(model1)
model2 <- weibreg(y ~ karno + diagtime + age + prior + celltype + trt, shape=0)
summary(model2)
model3 <- coxph(y ~ karno + diagtime + age + prior + celltype + trt)
summary(model3)

nobs <- nrow(y)
ftime <- veteran$time
ord <- order(ftime, decreasing=TRUE)
# ord <- 1:nobs
plot(0, type='n', xlim=c(0,1000), ylim=c(0,150), xlab='Follow-up time', ylab='Population')
segments(rep(0.0, nobs), 1:nobs, ftime[ord], 1:nobs, col='gray25')
cases <- veteran$status == 1
points((ftime[ord])[cases[ord]], (1:nobs)[cases[ord]], pch=20, col='red', cex=0.5)

m <- sum(cases) * 100
persons <- rep(1:nobs, rmultinom(1, m, ftime/sum(ftime)))
moments <- rep(NA, m)
for (i in 1:m) {
    moments[i] <- runif(1, min=0.0, max=(ftime)[persons[i]])
}
points(moments, match(persons, ord), pch=20, col='blue', cex=0.5)
legend('topright', legend=c('Case series', 'Base series'), col=c('red', 'blue'), pch=c(20, 20))

# Fit a logistic regression to the case-base sample:

predm <- cbind(karno, diagtime, age, prior == 10, celltype == 'squamous', celltype == 'smallcell', celltype == 'adeno', trt == 2)
head(predm)

d <- c(rep(0, m), rep(1, sum(cases)))
moments <- c(moments, ftime[cases])
persons <- c(persons, (1:nobs)[cases])
offset <- rep(log(sum(ftime)/m), length(d))

model4 <- glm(d ~ moments + karno[persons] + diagtime[persons] + age[persons] + prior[persons] +
             celltype[persons] + trt[persons] + offset(offset), family=binomial(link=logit))
summary(model4)
round(cbind(coef(model3), sqrt(diag(vcov(model3))),
            coef(model4)[3:length(coef(model4))], sqrt(diag(vcov(model4)))[3:length(coef(model4))]), 4)

lambda <- function(x, pred, beta) {
    return(as.numeric(exp(beta[1] + beta[2] * x + crossprod(beta[3:length(beta)], pred))))
}

surv <- rep(NA, nobs)
for (i in 1:nobs) {
    surv[i] <- exp(-integrate(lambda, lower=0, upper=90, pred=predm[i,], beta=coef(model4))$value)
}
mean(1.0 - surv)
sum(ftime <= 90)/nobs
plot(1.0 - surv, ftime, xlim=c(0,1))

# Fit a spline for time:
library(splines)
basis <- bs(moments)
model5 <- glm(d ~ basis + karno[persons] + diagtime[persons] + age[persons] + prior[persons] +
             celltype[persons] + trt[persons] + offset(offset), family=binomial(link=logit))
summary(model5)
round(cbind(coef(model3), sqrt(diag(vcov(model3))),
            coef(model5)[(2 + ncol(basis)):length(coef(model5))],
            sqrt(diag(vcov(model5)))[(2 + ncol(basis)):length(coef(model5))]), 4)

lambda <- function(x, pred, beta) {
    return(as.numeric(exp(beta[1] + as.numeric(predict(basis, x) %*% beta[2:(1 + ncol(basis))]) +
                          crossprod(beta[(2 + ncol(basis)):length(beta)], pred))))
}

spsurv <- rep(NA, nobs)
for (i in 1:nobs) {
    spsurv[i] <- exp(-integrate(lambda, lower=0, upper=90, pred=predm[i,], beta=coef(model5))$value)
}
mean(1.0 - spsurv)
sum(ftime <= 90)/nobs
plot(1.0 - spsurv, ftime, xlim=c(0,1))

plot(surv, spsurv)

# 4-5:

inpath <- '~/git_repositories/casebase/data/'
etimes <- read.csv(file.path(inpath, 'etimes.csv'), header=TRUE)
vtimes <- read.csv(file.path(inpath, 'vtimes.csv'), header=TRUE)
vtc <- vtimes[etimes[,1],2]
exposed <- (vtc < etimes[,2]) & (vtc >= etimes[,2] - 7.0)

popsize <- nrow(vtimes)
ncases <- nrow(etimes)
n.cases <- rep(1, ncases)
n.exposed.cases <- as.numeric(exposed)
n.unexposed.cases <- as.numeric(!exposed)

# Time-matched sampling:

sample.ctls = function(x) sample(1:popsize, x, FALSE)

ctl.case.ratio = 1
n.ctls = n.cases * ctl.case.ratio
denom.sample = lapply(n.ctls, sample.ctls)
nctls = sum(n.ctls)
n.exposed.ctls <- rep(NA, ncases)

set.seed(1)
cc <- 0
casestatus <- NULL
exposure <- NULL
matched <- NULL
for (i in 1:ncases) {
    cc <- cc + 1
    exstatus <- (vtimes[denom.sample[[i]],2] < etimes[i,2]) & (vtimes[denom.sample[[i]],2] >= etimes[i,2] - 7.0)
    n.exposed.ctls[i] = sum(exstatus)
    casestatus <- c(casestatus, rep(1, n.cases[i]), rep(0, n.ctls[i]))
    exposure <- c(exposure, rep(0, n.unexposed.cases[i]), rep(1, n.exposed.cases[i]), as.numeric(exstatus))
    matched <- c(matched, rep(cc, n.cases[i] + n.cases[i] * ctl.case.ratio))
}
# cbind(casestatus, exposure, matched)

n.unexposed.ctls <- n.ctls - n.exposed.ctls

model1 <- clogit(casestatus ~ exposure + strata(as.factor(matched)))
summary(model1)

maxage <- 140/7

# pdf(file.path(inpath, 'studybase.pdf'), width=8, height=6, paper='special')
postscript(file.path(inpath, 'studybase.eps'), width=8, height=6, paper='special', horizontal=FALSE)
op <- par(mar = c(4.25,4.25,1,1))
plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population')
rect(0,0,maxage,popsize, col="gray80",border=NA)
legend(0.025 * maxage, 0.975 * popsize, legend=c('Pop.-time'),
       pch=c(15), col=c('gray80'), pt.bg=c('gray80'),
       cex=1, bg='white', box.col='white')
par(op)
dev.off()

idx <- order(vtimes[,2])

# pdf(file.path(inpath, 'studybaseexposure.pdf'), width=8, height=6, paper='special')
postscript(file.path(inpath, 'studybaseexposure.eps'), width=8, height=6, paper='special', horizontal=FALSE)
op <- par(mar = c(4.25,4.25,1,1))
plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population')
rect(0,0,maxage,popsize, col="lightblue",border=NA)
segments(vtimes[idx,2]/7, (1:popsize), (vtimes[idx,2]+7)/7, (1:popsize), col='yellow2')
legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.-time','Exposed pop.-time'),
       pch=c(15,15), col=c('lightblue','yellow2'), pt.bg=c('lightblue','yellow2'),
       cex=1, bg='white', box.col='white')
par(op)
dev.off()

vidx <- order(vtc)
caseidx <- (1:popsize)[vtimes[idx,1] %in% etimes[,1]]

# pdf(file.path(inpath, 'caseseries.pdf'), width=8, height=6, paper='special')
postscript(file.path(inpath, 'caseseries.eps'), width=8, height=6, paper='special', horizontal=FALSE)
op <- par(mar = c(4.25,4.25,1,1))
plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population')
rect(0,0,maxage,popsize, col="lightblue",border=NA)
segments(vtimes[idx,2]/7, (1:popsize), (vtimes[idx,2]+7)/7, (1:popsize), col='yellow2')
points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75)
legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.-time','Exposed pop.-time','Case series'),
       pch=c(15,15,20), col=c('lightblue','yellow2','red'), pt.bg=c('lightblue','yellow2','red'),
       cex=1, bg='white', box.col='white')
par(op)
dev.off()

# pdf(file.path(inpath, 'baseseries.pdf'), width=8, height=6, paper='special')
postscript(file.path(inpath, 'baseseries.eps'), width=8, height=6, paper='special', horizontal=FALSE)
op <- par(mar = c(4.25,4.25,1,1))
plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population')
rect(0,0,maxage,popsize, col="lightblue",border=NA)
segments(vtimes[idx,2]/7, (1:popsize), (vtimes[idx,2]+7)/7, (1:popsize), col='yellow2')
points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75)
points(etimes[matched[casestatus==0],2]/7, vtimes[unlist(denom.sample),1], pch=20, col='blue', cex=0.75)
legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.-time','Exposed pop.-time','Case series','Base series'),
       pch=c(15,15,20,20), col=c('lightblue','yellow2','red','blue'), pt.bg=c('lightblue','yellow2','red','blue'),
       cex=1, bg='white', box.col='white')
par(op)
dev.off()

# Self-matched sampling:

minday <- 0
maxday <- 140

nbase <- 1
times <- rep(NA, 2 * (ncases + nbase * ncases))
casestatus <- rep(NA, 2 * (ncases + nbase * ncases))
exposure <- rep(NA, 2 * (ncases + nbase * ncases))
matched <- rep(NA, 2 * (ncases + nbase * ncases))

set.seed(1)
counter <- 0
casecounter <- 0
for (i in 1:ncases) {
    if (is.finite(etimes[i,2]) & etimes[i,2] > minday & etimes[i,2] <= maxday) {
        counter <- counter + 1
        casecounter <- casecounter + 1
        times[counter] <- etimes[i,2]
        casestatus[counter] <- 1
        exposure[counter] <- (etimes[i,2] - vtc[i]) > 0.0 & (etimes[i,2] - vtc[i]) < 7.0
        matched[counter] <- casecounter
        # n <- rpois(1, nbase)
        n <- nbase
        for (j in 1:n) {
            counter <- counter + 1
            u <- runif(1, 0.0, nbase)
            times[counter] <- runif(1, min=minday, max=maxday)
            casestatus[counter] <- 0
            exposure[counter] <- (times[counter] - vtc[i]) > 0.0 & (times[counter] - vtc[i]) < 7.0
            matched[counter] <- casecounter
        }
    }
    # if (i %% 1000 == 0)
    #     print(i)
}
# counter
# casecounter
times <- times[1:counter]
casestatus <- casestatus[1:counter]
exposure <- exposure[1:counter]
matched <- matched[1:counter]

times2 <- times^2
model2 <- clogit(casestatus ~ exposure + times + times2 + strata(as.factor(matched)))
summary(model2)

model3 <- clogit(casestatus ~ exposure + bs(times) + strata(as.factor(matched)))
summary(model3)

# pdf(file.path(inpath, 'caseseries2.pdf'), width=8, height=6, paper='special')
postscript(file.path(inpath, 'caseseries2.eps'), width=8, height=6, paper='special', horizontal=FALSE)
op <- par(mar = c(4.25,4.25,1,1))
plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population')
rect(0,0,maxage,popsize, col="lightblue",border=NA)
segments(vtimes[idx,2]/7, (1:popsize), (vtimes[idx,2]+7)/7, (1:popsize), col='yellow2')
points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75)
legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.n-time','Exposed pop.-time','Case series'),
       pch=c(15,15,20), col=c('lightblue','yellow2','red'), pt.bg=c('lightblue','yellow2','red'),
       cex=1, bg='white', box.col='white')
par(op)
dev.off()

# pdf(file.path(inpath, 'selfmatched.pdf'), width=8, height=6, paper='special')
postscript(file.path(inpath, 'selfmatched.eps'), width=8, height=6, paper='special', horizontal=FALSE)
op <- par(mar = c(4.25,4.25,1,1))
plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population')
rect(0,0,maxage,popsize, col="white",border=NA)
segments(rep(minday/7, ncases), caseidx, rep(maxday/7, ncases), caseidx, col='lightblue')
segments(vtc[vidx]/7, caseidx, (vtc[vidx]+7)/7, caseidx, col='yellow2')
points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75)
legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.-time','Exposed pop.-time','Case series'),
       pch=c(15,15,20), col=c('lightblue','yellow2','red'), pt.bg=c('lightblue','yellow2','red'),
       cex=1, bg='white', box.col='white')
par(op)
dev.off()

# pdf(file.path(inpath, 'baseseries2.pdf'), width=8, height=6, paper='special')
postscript(file.path(inpath, 'baseseries2.eps'), width=8, height=6, paper='special', horizontal=FALSE)
op <- par(mar = c(4.25,4.25,1,1))
plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population')
rect(0,0,maxage,popsize, col="white",border=NA)
segments(rep(minday/7, ncases), caseidx, rep(maxday/7, ncases), caseidx, col='lightblue')
segments(vtc[vidx]/7, caseidx, (vtc[vidx]+7)/7, caseidx, col='yellow2')
points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75)
points(times[casestatus==0]/7, caseidx[matched[casestatus==0]], pch=20, col='blue', cex=0.75)
legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.-time','Exposed pop.-time','Case series','Base series'),
       pch=c(15,15,20,20), col=c('lightblue','yellow2','red','blue'), pt.bg=c('lightblue','yellow2','red','blue'),
       cex=1, bg='white', box.col='white')
par(op)
dev.off()

# pdf(file.path(inpath, 'baseseries3.pdf'), width=8, height=6, paper='special')
postscript(file.path(inpath, 'baseseries3.eps'), width=8, height=6, paper='special', horizontal=FALSE)
op <- par(mar = c(4.25,4.25,1,1))
plot(1, 1, xlim=c(0,maxage), ylim=c(0,1000), type='n', xlab='Age (weeks)', ylab='Population')
rect(0,0,maxage,popsize, col="white",border=NA)
segments(rep(minday/7, ncases), caseidx, rep(maxday/7, ncases), caseidx, col='lightblue')
segments(vtc[vidx]/7, caseidx, (vtc[vidx]+7)/7, caseidx, col='yellow2')
points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75)
points(times[casestatus==0]/7, caseidx[matched[casestatus==0]], pch=20, col='blue', cex=0.75)
par(op)
dev.off()


