#R code, updated 2014-03-26
library(magrittr)
# library(casebase)
library(data.table)

devtools::load_all()
dev.off()
rm(list=ls())
# setwd("~/Dropbox/work/compete/data")

PATIENT=read.table("data-raw/olli/patient.csv", sep=",", head=T, stringsAsFactors = FALSE)
pt=PATIENT[,c("pid", "rndgroup", "age", "death_days",
"deathcutoff", "finaldeathLC","fup_days","evpdeath","candx_days")]
dim(pt); length(pt$pid)
str(pt)
head(pt)

pt$fuptime=as.numeric(as.character(pt$fup_days))
summary(pt$fuptime/365)

pt$deathtime=as.numeric(as.character(pt$death_days))
summary(pt$deathtime/365)

pt$cancertime=as.numeric(as.character(pt$candx_days))
pt$cancertime=ifelse(pt$cancertime==0, 1, pt$cancertime)
summary(pt$cancertime/365)

head(pt)
pt$pid %>% unique %>% length()
nrow(pt)

head(pt)

DT_pt <- as.data.table(pt)
setkey(DT_pt, "pid")

DT_pt[, event:= 1 * (!is.na(deathtime) & finaldeathLC %in% c('1')) +
       2 * (!is.na(deathtime) & finaldeathLC %in% c('0','M')), by = pid]


DT_pt[, table(event, rndgroup)]
DT_pt[, rndgroup:=factor(rndgroup, levels = 1:2, labels = c("X-ray arm","CT arm"))]


popTimeObj <- popTime(DT_pt, time = "fuptime", event = "event", exposure = "rndgroup")
p5 <- plot(popTimeObj)
str(p5)


popTimeObj$data$ycoord %>% max
p5  + scale_y_continuous(breaks = seq(0, max(popTimeObj$data$ycoord), 5000))


devtools::load_all()

DT <- read.csv(system.file("extdata", "bmtcrr.csv", package = "casebase"))
popTimeData <- popTime(data = DT, time = "ftime", exposure = "D")
p <- plot(popTimeData)
p + scale_y_continuous(breaks = seq(0, max(popTimeData$data$ycoord), 10))

popTimeData <- popTime(data = veteran)

pt0=subset(pt,pt$rndgroup==2);dim(pt0)
n0 <- nrow(pt0)
pt1=subset(pt,pt$rndgroup==1);dim(pt1)
n1 <- nrow(pt1)
sum(pt0$fuptime/365.25)
sum(pt1$fuptime/365.25)

# Lung cancer mortality outcome:

table(!is.na(pt$deathtime), pt$finaldeathLC)
table(pt$finaldeathLC, pt$evpdeath)

#pt0_sort=pt0[order(pt0$fuptime,decreasing=TRUE),]
#pt1_sort=pt1[order(pt1$fuptime,decreasing=TRUE),]

evtype0 = 1 * (!is.na(pt0$deathtime) & pt0$finaldeathLC %in% c('1')) +
       2 * (!is.na(pt0$deathtime) & pt0$finaldeathLC %in% c('0','M'))
table(evtype0)
evtype1 = 1 * (!is.na(pt1$deathtime) & pt1$finaldeathLC %in% c('1')) +
       2 * (!is.na(pt1$deathtime) & pt1$finaldeathLC %in% c('0','M'))
table(evtype1)
when.event0=pt0$deathtime[evtype0==1]
when.event1=pt1$deathtime[evtype1==1]

evtime0 <- pt0$fuptime
evtime1 <- pt1$fuptime

ncc <- length(unique(c(evtype0,evtype1)))

Day0=1:(max(evtime0)+1)
N0=matrix(NA, ncc, max(evtime0)+1)
for (d in 1:(max(evtime0)+1)) {
    for (c in 1:ncc) {
	    N0[c,d]= sum(evtime0[evtype0==(c-1)]<d)
    }
}

Day1=1:(max(evtime1)+1)
N1=matrix(NA, ncc, max(evtime1)+1)
for (d in 1:(max(evtime1)+1)) {
    for (c in 1:ncc) {
	    N1[c,d]= sum(evtime1[evtype1==(c-1)]<d)
    }
}

y.d0 = N0[2,when.event0] + (n0 - (colSums(N0[,when.event0]))) * runif(sum(evtype0==1))
y.d1 = N1[2,when.event1] + (n1 - (colSums(N1[,when.event1]))) * runif(sum(evtype1==1))

pdf("~/Dropbox/work/compete/NLST_cancer_mortality.pdf", width=9.5, height=9, paper="special")
op <- par(mar=c(6, 5, 2, 2) + 0.1)
N=nrow(pt);max.yr=max(pt$fuptime/365.25)
plot(1,1, xlab="Follow-up Year",ylab="",ylim=c(0,N),xlim=c(0,max.yr),type="n",yaxt="n",cex.lab=1.2)
mtext("Population", side = 2, cex=1.2,line=4)
rect(0,0,max.yr,n0,col="grey90",border=NA)
rect(0,n0,max.yr,N,col="yellow2",border=NA)

polygon(c(Day0/365.25,max.yr,max.yr,0),c(n0-(N0[1,]+N0[3,]),n0-max(N0[1,]+N0[3,]),n0,n0),col='cyan',border=NA)
polygon(c(Day0/365.25,max.yr,max.yr,0),c(n0-N0[1,],n0-max(N0[1,]),n0,n0),col='black',border=NA)
polygon(c(Day0/365.25,max.yr,max.yr,0),c(N0[2,],max(N0[2,]),0,0),col='red',border=NA)
points(when.event0/365.25,y.d0, col='red',pch=16,cex=.5)

polygon(c(Day1/365.25,max.yr,max.yr,0),n0+c(n1-(N1[1,]+N1[3,]),n1-max(N1[1,]+N1[3,]),n1,n1),col='cyan',border=NA)
polygon(c(Day1/365.25,max.yr,max.yr,0),n0+c(n1-N1[1,],n1-max(N1[1,]),n1,n1),col='black',border=NA)
polygon(c(Day1/365.25,max.yr,max.yr,0),n0+c(N1[2,],max(N1[2,]),0,0),col='red',border=NA)
points(when.event1/365.25, n0+y.d1, col='red',pch=16,cex=.5)

abline(h=n0,col="white")
text(7.5,n1+n0/2,"CT arm",col='white')
text(7.5,n0/2,"X-ray arm",col='white')

legend(-1, -5000, legend=c("Lung cancer mortality","Other mortality","Censoring"),pch=c(15,15,15),col=c('red','cyan','black'), bty="n", xpd=TRUE)
axis(2,at=seq(0,n0,by=5000),labels=seq(0,n0,by=5000),las=2,hadj=0.85)
axis(2,at=n0+seq(0,n1,by=5000),labels=seq(0,n0,by=5000),las=2,hadj=0.85)
dev.off()

# Cancer incidence outcome:

table(!is.na(pt0$deathtime), !is.na(pt0$cancertime))
table(pt0$deathtime < pt0$cancertime, useNA='ifany')

table(!is.na(pt1$deathtime), !is.na(pt1$cancertime))
table(pt1$deathtime < pt1$cancertime, useNA='ifany')

evtype0 = 1 * (!is.na(pt0$cancertime)) + 2 * (!is.na(pt0$deathtime) & is.na(pt0$cancertime))
table(evtype0)
evtime0 <- pmin(pt0$deathtime, pt0$cancertime, pt0$fuptime, na.rm=TRUE)

evtype1 = 1 * (!is.na(pt1$cancertime)) + 2 * (!is.na(pt1$deathtime) & is.na(pt1$cancertime))
table(evtype1)
evtime1 <- pmin(pt1$deathtime, pt1$cancertime, pt1$fuptime, na.rm=TRUE)

when.event0=evtime0[evtype0==1]
when.event1=evtime1[evtype1==1]

ncc <- length(unique(c(evtype0,evtype1)))

Day0=1:(max(evtime0)+1)
N0=matrix(NA, ncc, max(evtime0)+1)
for (d in 1:(max(evtime0)+1)) {
    for (c in 1:ncc) {
	    N0[c,d]= sum(evtime0[evtype0==(c-1)]<d)
    }
}

Day1=1:(max(evtime1)+1)
N1=matrix(NA, ncc, max(evtime1)+1)
for (d in 1:(max(evtime1)+1)) {
    for (c in 1:ncc) {
	    N1[c,d]= sum(evtime1[evtype1==(c-1)]<d)
    }
}

y.d0 = N0[2,when.event0] + (n0 - (colSums(N0[,when.event0]))) * runif(sum(evtype0==1))
y.d1 = N1[2,when.event1] + (n1 - (colSums(N1[,when.event1]))) * runif(sum(evtype1==1))

pdf("~/Dropbox/work/compete/NLST_cancer_incidence.pdf", width=9.5, height=9, paper="special")
op <- par(mar=c(6, 5, 2, 2) + 0.1)
N=nrow(pt);max.yr=max(pt$fuptime/365.25)
plot(1,1, xlab="Follow-up Year",ylab="",ylim=c(0,N),xlim=c(0,max.yr),type="n",yaxt="n",cex.lab=1.2)
mtext("Population", side = 2, cex=1.2,line=4)
rect(0,0,max.yr,n0,col="grey90",border=NA)
rect(0,n0,max.yr,N,col="yellow2",border=NA)

polygon(c(Day0/365.25,max.yr,max.yr,0),c(n0-(N0[1,]+N0[3,]),n0-max(N0[1,]+N0[3,]),n0,n0),col='cyan',border=NA)
polygon(c(Day0/365.25,max.yr,max.yr,0),c(n0-N0[1,],n0-max(N0[1,]),n0,n0),col='black',border=NA)
polygon(c(Day0/365.25,max.yr,max.yr,0),c(N0[2,],max(N0[2,]),0,0),col='red',border=NA)
points(when.event0/365.25,y.d0, col='red',pch=16,cex=.5)

polygon(c(Day1/365.25,max.yr,max.yr,0),n0+c(n1-(N1[1,]+N1[3,]),n1-max(N1[1,]+N1[3,]),n1,n1),col='cyan',border=NA)
polygon(c(Day1/365.25,max.yr,max.yr,0),n0+c(n1-N1[1,],n1-max(N1[1,]),n1,n1),col='black',border=NA)
polygon(c(Day1/365.25,max.yr,max.yr,0),n0+c(N1[2,],max(N1[2,]),0,0),col='red',border=NA)
points(when.event1/365.25, n0+y.d1, col='red',pch=16,cex=.5)

abline(h=n0,col="white")
text(7.5,n1+n0/2,"CT arm",col='white')
text(7.5,n0/2,"X-ray arm",col='white')

legend(-1, -5000, legend=c("Lung cancer incidence","Mortality","Censoring"),pch=c(15,15,15),col=c('red','cyan','black'), bty="n", xpd=TRUE)
axis(2,at=seq(0,n0,by=5000),labels=seq(0,n0,by=5000),las=2,hadj=0.85)
axis(2,at=n0+seq(0,n1,by=5000),labels=seq(0,n0,by=5000),las=2,hadj=0.85)
dev.off()



