# jh july 27, 2016


# setwd("/Users/jameshanley/Dropbox/work/osm/profileArticle")

# ds=read.csv("ERSPCindividualData.csv")
rm(list=ls())
# ERSPC <- read.csv("data-raw/hanley/ERSPCindividualData.csv")
# devtools::use_data(ERSPC, overwrite = TRUE)
devtools::load_all()
data("ERSPC")

# load(file = "data/ERSPC.rda")

# str(ds); head(ds)

library(survival)
library(casebase)
library(data.table)

# DT_ds <- as.data.table(ds)
# DT_ds[, ScrArm:=factor(ScrArm, levels = 0:1, labels = c("No-Screening Arm","Screening Arm"))]
# popTimeData <- popTime(DT_ds, event = "DeadOfPrCa", exposure = "ScrArm")
# plot(popTimeData)


# cumulative incidence plots
ds <- ERSPC
KM <- survfit(Surv(Follow.Up.Time,DeadOfPrCa) ~ ScrArm, data = ds)
str(KM)

par(mfrow=c(1,1),mar = c(5,5,0.1,0.1))
plot(KM$time[    1: 1501], 1-KM$surv[   1:1501], type="s", col="red" ,
  ylab = "Risk", xlab="Years since Randomization")
lines(KM$time[1502: 2923], 1-KM$surv[1502: 2923], type="s", col="green" )


#  PopulationTime plots


par(mfrow=c(1,1),mar = c(0.01,0.01,0.1,0.1))

plot(c(-0.5,15.75),c(-93000,80000), col="white" )
set.seed(7654321)

OFF = 2000


for(i in 0:1) {
	t=seq(0.01,14.9,0.01)
	S = function(x) sum(ds$Follow.Up.Time[ds$ScrArm==i] >= x)
	n = unlist(lapply(t,"S"))
	if(i==1) yy =  c(0,n,0) + OFF
	if(i==0) yy =  c(0,-n,0) - OFF
	polygon(c(0,t,14.9),yy,col="grey80",border=NA)

    t.d = ds$Follow.Up.Time[ds$ScrArm==i & ds$DeadOfPrCa==1]

	for( j in 1:length(t.d) ) {
		time.index =  ceiling(t.d[j]/0.01)
		nn   = n[ time.index ]
		if(i==1) h = runif(1,0.01*nn,0.99*nn)  + OFF
		if(i==0) h = runif(1,-0.99*nn,-0.01*nn) - OFF
		points(t.d[j],h, pch=19,cex=0.25,col="red")
	}
}

for (t in 1:15) text(t,0,toString(t), cex=0.75)
text(15.25,0,"Year", cex=0.75,adj=c(0,0.5))

for (n in seq(0,90000,10000)) {
	if(n> 0 & n < 80000) text(-0.1,n+OFF,format(n,big.mark=","), cex=0.75,adj=c(1,0.5))
	if(n> 0) text(-0.1,-n-OFF,format(n,big.mark=","), cex=0.75,adj=c(1,0.5))
	segments(-0.05,  n+OFF, 0, n+OFF , lwd=0.5)
	segments(-0.05, -n-OFF, 0, -n-OFF, lwd=0.5 )

}
text(4, 70000+OFF,"Screening Arm of ERSPC", cex=1,adj=c(0,0.5))
text(4,-85000-OFF,"No-Screening Arm", cex=1,adj=c(0,0.5))

text(-0.75,78000+OFF,"Number of
Men being Followed", cex=1,adj=c(0,0.5))
h = 50000+OFF
points(9.5,h, pch=19,cex=0.25,col="red")
text(9.6,h,"Death from Prostate Cancer", adj=c(0,0.5))



