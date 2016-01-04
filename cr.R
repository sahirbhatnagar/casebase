set.seed(1)
#setwd('~/Dropbox/work/CHL5209H_2015/slides')
#getwd()

N <- 10000
a1 <- 1
a2 <- 500
b1 <- 1
b2 <- 5000

x <- rbinom(N, 1, 0.5)
z <- rbinom(N, 1, 0.5)
table(x, z)

mi <- rweibull(N, a1, a2 * exp(2 * z + 2 * x)^(-1/a1))
st <- rweibull(N, b1, b2 * exp(2 * z + 2 * x)^(-1/b1))

yrs <- 10
mie <- mi < yrs
table(mie)
ste <- st < yrs
table(ste)
mi[!mie] <- yrs
st[!ste] <- yrs

stc <- st < mi
table(stc)
mic <- mi < st
table(mic)
t <- pmin(st, mi)

ycoord <- rep(NA, N)
for (i in 1:N) {
    ycoord[i] <- runif(1, 0, sum(mi >= st[i]))
}
plot(st, ycoord)

pdf(file.path(getwd(), 'studybase.pdf'), height=7, width=7, paper='special')
idx <- order(!mic, ifelse(mic, t, 0), decreasing=TRUE)
plot(1,1,type='n', xlim=c(0,yrs), ylim=c(1,N), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, N), (1:N), rep(yrs, N), (1:N), col='lightblue')
dev.off()

pdf(file.path(getwd(), 'studybase2.pdf'), height=7, width=7, paper='special')
idx <- order(!mic, ifelse(mic, t, 0), decreasing=TRUE)
plot(1,1,type='n', xlim=c(0,yrs), ylim=c(1,N), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, N), (1:N), rep(yrs, N), (1:N), col='lightblue')
op <- par(lwd=3)
legend(0.1, 100, legend=c('Person-time'), lty=c('solid'), col=c('lightblue'), bg='white', xjust=0, yjust=0, box.lwd=1)
par(op)
dev.off()

pdf(file.path(getwd(), 'atrisk.pdf'), height=7, width=7, paper='special')
idx <- order(!mic, ifelse(mic, t, 0), decreasing=TRUE)
plot(1,1,type='n', xlim=c(0,yrs), ylim=c(1,N), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, N), (1:N), rep(yrs, N), (1:N), col='gray60')
segments(rep(0.0, N), (1:N), t[idx], (1:N), col='lightblue')
points((t[idx])[mic[idx]], (1:N)[mic[idx]], pch=20, col='red', cex=0.5)
text(7.5, 9000, expression(N[MI](t)==1), cex=1.5)
op <- par(lwd=3)
legend(0.1, 100, legend=c('Incident MI event','At-risk person-time'), pch=c(20, NA), lty=c(NA, 'solid'), col=c('red','lightblue'), bg='white',  pt.cex=c(1, NA), xjust=0, yjust=0, box.lwd=1)
par(op)
# points((t[idx])[stc[idx]], (1:N)[stc[idx]], pch=20, col='blue')
dev.off()

pdf(file.path(getwd(), 'events.pdf'), height=7, width=7, paper='special')
idx <- order(!mic, ifelse(mic, t, 0), decreasing=TRUE)
plot(1,1,type='n', xlim=c(0,yrs), ylim=c(1,N), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, N), (1:N), rep(yrs, N), (1:N), col='gray60')
segments(rep(0.0, N), (1:N), t[idx], (1:N), col='lightblue')
points((t[idx])[mic[idx]], (1:N)[mic[idx]], pch=20, col='red', cex=0.5)
points((t[idx])[stc[idx]], (ycoord[idx])[stc[idx]], pch=20, col='blue', cex=0.5)
text(7.5, 9000, expression(N[MI](t)==1), cex=1.5)
op <- par(lwd=3)
legend(0.1, 100, legend=c('Incident stroke event','Incident MI event','At-risk person-time'), pch=c(20, 20, NA), lty=c(NA, NA, 'solid'), col=c('blue', 'red','lightblue'), bg='white',  pt.cex=c(1, 1, NA), xjust=0, yjust=0, box.lwd=1)
par(op)
dev.off()

pdf(file.path(getwd(), 'byz.pdf'), height=7, width=7, paper='special')
idx <- order(z, !mic, ifelse(mic, t, 0), decreasing=TRUE)
plot(1,1,type='n', xlim=c(0,yrs), ylim=c(1,N), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, N), (1:N), rep(yrs, N), (1:N), col='gray60')
segments(rep(0.0, N), (1:N), t[idx], (1:N), col='lightblue')
segments(rep(0.0, N), (1:N)[z[idx]==1], (t[idx])[z[idx]==1], (1:N)[z[idx]==1], col='yellow2')
points((t[idx])[mic[idx]], (1:N)[mic[idx]], pch=20, col='red', cex=0.5)
 # text(8.5, 9250, expression(N[MI](t)==1), cex=1.5)
text(8.5, 4250, expression(N[MI](t)==1), cex=1.5)
text(2, 2500, expression(Z(t)==1), cex=1.5)
text(2, 7500, expression(Z(t)==0), cex=1.5)
op <- par(lwd=3)
legend(0.1, 100, legend=c('Non-smoking person-time','Smoking person-time'), lty=c('solid','solid'),
col=c('lightblue','yellow2'), bg='white', xjust=0, yjust=0, box.lwd=1)
par(op)
dev.off()

ycoord <- rep(NA, N)
for (i in 1:N) {
    if (z[i] == 1)
        ycoord[i] <- runif(1, 0, sum(mi >= st[i] & z == z[i]))
    else if (z[i] == 0)
        ycoord[i] <- runif(1, sum(z == 1), sum(z == 1) + sum(mi >= st[i] & z == z[i]))
}
plot(st, ycoord)

pdf(file.path(getwd(), 'byzevents.pdf'), height=7, width=7, paper='special')
idx <- order(z, !mic, ifelse(mic, t, 0), decreasing=TRUE)
plot(1,1,type='n', xlim=c(0,yrs), ylim=c(1,N), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, N), (1:N), rep(yrs, N), (1:N), col='gray60')
segments(rep(0.0, N), (1:N), t[idx], (1:N), col='lightblue')
segments(rep(0.0, N), (1:N)[z[idx]==1], (t[idx])[z[idx]==1], (1:N)[z[idx]==1], col='yellow2')
points((t[idx])[mic[idx]], (1:N)[mic[idx]], pch=20, col='red', cex=0.5)
points((t[idx])[stc[idx]], (ycoord[idx])[stc[idx]], pch=20, col='blue', cex=0.5)
 # text(8.5, 9250, expression(N[MI](t)==1), cex=1.5)
text(8.5, 4250, expression(N[MI](t)==1), cex=1.5)
dev.off()

pdf(file.path(getwd(), 'byzrates.pdf'), height=7, width=7, paper='special')
idx <- order(z, !mic, ifelse(mic, t, 0), decreasing=TRUE)
plot(1,1,type='n', xlim=c(0,yrs), ylim=c(1,N), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, N), (1:N), rep(yrs, N), (1:N), col='gray60')
segments(rep(0.0, N), (1:N), t[idx], (1:N), col='lightblue')
segments(rep(0.0, N), (1:N)[z[idx]==1], (t[idx])[z[idx]==1], (1:N)[z[idx]==1], col='yellow2')
points((t[idx])[mic[idx]], (1:N)[mic[idx]], pch=20, col='red', cex=0.5)
 # text(8.5, 9250, expression(N[MI](t)==1), cex=1.5)
text(8.5, 4250, expression(N[MI](t)==1), cex=1.5)
text(5, 2000, expression(hat(lambda)[S](paste(t, '|', Z(t)==1)) == frac(184, 38311) %~~% frac(4.8, 1000)), cex=1.5)
text(5, 7000, expression(hat(lambda)[S](paste(t, '|', Z(t)==0)) == frac(54, 47614) %~~% frac(1.1, 1000)), cex=1.5)
dev.off()

pdf(file.path(getwd(), 'byzx.pdf'), height=7, width=7, paper='special')
idx <- order(z, x, !mic, ifelse(mic, t, 0), decreasing=TRUE)
plot(1,1,type='n', xlim=c(0,yrs), ylim=c(1,N), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, N), (1:N), rep(yrs, N), (1:N), col='gray60')
segments(rep(0.0, N), (1:N), t[idx], (1:N), col='lightblue')
segments(rep(0.0, N), (1:N)[z[idx]==1], (t[idx])[z[idx]==1], (1:N)[z[idx]==1], col='yellow2')
points((t[idx])[mic[idx]], (1:N)[mic[idx]], pch=20, col='red', cex=0.5)
 # text(8.5, 9250, expression(N[MI](t)==1), cex=1.5)
text(2, 500, expression(paste(Z(t)==1,', ',X(t)==1)), cex=1.25)
text(2, 3000, expression(paste(Z(t)==1,', ',X(t)==0)), cex=1.25)
text(2, 5500, expression(paste(Z(t)==0,', ',X(t)==1)), cex=1.25)
text(2, 8000, expression(paste(Z(t)==0,', ',X(t)==0)), cex=1.25)
text(8.5, 1900, expression(N[MI](t)==1), cex=1.25)
dev.off()

pdf(file.path(getwd(), 'byzxrates.pdf'), height=7, width=7, paper='special')
idx <- order(z, x, !mic, ifelse(mic, t, 0), decreasing=TRUE)
plot(1,1,type='n', xlim=c(0,yrs), ylim=c(1,N), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, N), (1:N), rep(yrs, N), (1:N), col='gray60')
segments(rep(0.0, N), (1:N), t[idx], (1:N), col='lightblue')
segments(rep(0.0, N), (1:N)[z[idx]==1], (t[idx])[z[idx]==1], (1:N)[z[idx]==1], col='yellow2')
points((t[idx])[mic[idx]], (1:N)[mic[idx]], pch=20, col='red', cex=0.5)
 # text(8.5, 9250, expression(N[MI](t)==1), cex=1.5)
text(4, 600, expression(hat(lambda)[S](paste(t, '|', Z(t)==1, ', ', X(t)==1)) == frac(154, 14878) %~~% frac(10.2, 1000)), cex=1.25)
text(4, 3100, expression(hat(lambda)[S](paste(t, '|', Z(t)==1, ', ', X(t)==0)) == frac(32, 23433) %~~% frac(1.4, 1000)), cex=1.25)
text(4, 5600, expression(hat(lambda)[S](paste(t, '|', Z(t)==0, ', ', X(t)==1)) == frac(48, 22647) %~~% frac(2.1, 1000)), cex=1.25)
text(4, 8100, expression(hat(lambda)[S](paste(t, '|', Z(t)==0, ', ', X(t)==0)) == frac(5, 24967) %~~% frac(0.2, 1000)), cex=1.25)
text(8.5, 1900, expression(N[MI](t)==1), cex=1.25)
dev.off()

ycoord <- rep(NA, N)
for (i in 1:N) {
    if (z[i] == 1 & x[i] == 1)
        ycoord[i] <- runif(1, 0, sum(mi >= st[i] & z == z[i] & x == x[i]))
    else if (z[i] == 1 & x[i] == 0)
        ycoord[i] <- runif(1, sum(z == 1 & x == 1), sum(z == 1 & x == 1) + sum(mi >= st[i] & z == z[i] & x == x[i]))
    else if (z[i] == 0 & x[i] == 1)
        ycoord[i] <- runif(1, sum(z == 1), sum(z == 1) + sum(mi >= st[i] & z == z[i] & x == x[i]))
    else if (z[i] == 0 & x[i] == 0)
        ycoord[i] <- runif(1, sum(!(z == 0 & x == 0)), sum(!(z == 0 & x == 0)) + sum(mi >= st[i] & z == z[i] & x == x[i]))
}
plot(st, ycoord)

pdf(file.path(getwd(), 'byzxevents.pdf'), height=7, width=7, paper='special')
idx <- order(z, x, !mic, ifelse(mic, t, 0), decreasing=TRUE)
plot(1,1,type='n', xlim=c(0,yrs), ylim=c(1,N), xlab='Follow-up years', ylab='Population')
segments(rep(0.0, N), (1:N), rep(yrs, N), (1:N), col='gray60')
segments(rep(0.0, N), (1:N), t[idx], (1:N), col='lightblue')
segments(rep(0.0, N), (1:N)[z[idx]==1], (t[idx])[z[idx]==1], (1:N)[z[idx]==1], col='yellow2')
points((t[idx])[mic[idx]], (1:N)[mic[idx]], pch=20, col='red', cex=0.5)
points((t[idx])[stc[idx]], (ycoord[idx])[stc[idx]], pch=20, col='blue', cex=0.5)
 # text(8.5, 9250, expression(N[MI](t)==1), cex=1.5)
text(8.5, 1900, expression(N[MI](t)==1), cex=1.25)
dev.off()


1000 * sum(stc[z==1 & x==1])/sum(t[z==1 & x==1])
1000 * sum(stc[z==1 & x==0])/sum(t[z==1 & x==0])
1000 * sum(stc[z==0 & x==1])/sum(t[z==0 & x==1])
1000 * sum(stc[z==0 & x==0])/sum(t[z==0 & x==0])

sum(ste)/sum(st)
(sum(ste[z==1])/sum(st[z==1]))/(sum(ste[z==0])/sum(st[z==0]))

sum(stc)/sum(t)
(sum(stc[z==1])/sum(t[z==1]))/(sum(stc[z==0])/sum(t[z==0]))



