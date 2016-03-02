---
title: "Fitting smooth-in-time parametric hazard functions"
author: "Maxime Turgeon"
date: "2016-02-29"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: reference.bib
---

## Methodological details

Case-base sampling was proposed by Hanley & Miettinen [-@hanley2009fitting] as a way to fit smooth-in-time parametric hazard functions via logistic regression. The main idea, which was first proposed by Mantel [-@mantel1973synthetic] and then later developped by Efron [-@efron1977efficiency], is to sample person-moments, i.e. discrete time points along an subject's follow-up time, in order to construct a base series against which the case series can be compared. 

This approach allows the explicit inclusion of the time variable into the model, which enables the user to fit a wide class of parametric hazard functions. For example, including time linearly recovers the Gompertz hazard, whereas including time *logarithmically* recovers the Weibull hazard; not including time at all corresponds to the exponential hazard.

The theoretical properties of this approach have been studied in Saarela & Arjas [-@saarela2015non] and Saarela [-@saarela2015case].

## First example

The first example we discuss uses the well-known ```veteran``` dataset, which is part of the ```survival``` package. As we can see below, there is almost no censoring, and therefore we can get a good visual representation of the survival function:


```r
library(survival)
```

```
## 
## Attaching package: 'survival'
```

```
## The following object is masked _by_ '.GlobalEnv':
## 
##     veteran
```

```r
data(veteran)
table(veteran$status)
```

```
## 
##   0   1 
##   9 128
```

```r
evtimes <- veteran$time[veteran$status == 1]
hist(evtimes, nclass=30, main='', xlab='Survival time (days)', col='gray90', probability=TRUE)
tgrid <- seq(0, 1000, by=10)
lines(tgrid, dexp(tgrid, rate=1.0/mean(evtimes)), 
      lwd=2, lty=2, col='red')
```

<img src="smoothHazard_files/figure-html/unnamed-chunk-1-1.png" title="" alt="" style="display: block; margin: auto;" />

As we can see, the empirical survival function ressembles an exponential distribution.

We will first try to estimate the hazard function parametrically using some well-known regression routines. But first, we will reformat the data slightly.


```r
veteran$prior <- factor(veteran$prior, levels = c(0, 10))
veteran$celltype <- factor(veteran$celltype, 
                           levels = c('large', 'squamous', 'smallcell', 'adeno'))
veteran$trt <- factor(veteran$trt, levels = c(1, 2))
```

Using the ```eha``` package, we can fit a Weibull form, with different values of the shape parameter. For $shape = 1$, we get an exponential distribution:


```r
library(eha)
```

```
## 
## Attaching package: 'eha'
```

```
## The following objects are masked from 'package:VGAM':
## 
##     dgompertz, dmakeham, pgompertz, pmakeham, qgompertz, qmakeham,
##     rgompertz, rmakeham
```

```r
y <- with(veteran, Surv(time, status))

model1 <- weibreg(y ~ karno + diagtime + age + prior + celltype + trt, 
                  data = veteran, shape = 1)
summary(model1)
```

```
## Call:
## weibreg(formula = y ~ karno + diagtime + age + prior + celltype + 
##     trt, data = veteran, shape = 1)
## 
## Covariate           Mean       Coef Exp(Coef)  se(Coef)    Wald p
## karno              68.419    -0.031     0.970     0.005     0.000 
## diagtime            8.139     0.000     1.000     0.009     0.974 
## age                57.379    -0.006     0.994     0.009     0.505 
## prior 
##                0    0.653     0         1           (reference)
##               10    0.347     0.049     1.051     0.227     0.827 
## celltype 
##            large    0.269     0         1           (reference)
##         squamous    0.421    -0.377     0.686     0.273     0.166 
##        smallcell    0.206     0.443     1.557     0.261     0.090 
##            adeno    0.104     0.736     2.087     0.294     0.012 
## trt 
##                1    0.477     0         1           (reference)
##                2    0.523     0.220     1.246     0.199     0.269 
## 
## log(scale)                    2.811    16.633     0.713     0.000 
## 
##  Shape is fixed at  1 
## 
## Events                    128 
## Total time at risk         16663 
## Max. log. likelihood      -716.16 
## LR test statistic         70.1 
## Degrees of freedom        8 
## Overall p-value           4.64229e-12
```

If we take $shape = 0$, the shape parameter is estimated along with the regression coefficients:


```r
model2 <- weibreg(y ~ karno + diagtime + age + prior + celltype + trt, 
                  data = veteran, shape = 0)
summary(model2)
```

```
## Call:
## weibreg(formula = y ~ karno + diagtime + age + prior + celltype + 
##     trt, data = veteran, shape = 0)
## 
## Covariate           Mean       Coef Exp(Coef)  se(Coef)    Wald p
## karno              68.419    -0.032     0.968     0.005     0.000 
## diagtime            8.139     0.001     1.001     0.009     0.955 
## age                57.379    -0.007     0.993     0.009     0.476 
## prior 
##                0    0.653     0         1           (reference)
##               10    0.347     0.047     1.048     0.229     0.836 
## celltype 
##            large    0.269     0         1           (reference)
##         squamous    0.421    -0.428     0.651     0.278     0.123 
##        smallcell    0.206     0.462     1.587     0.262     0.078 
##            adeno    0.104     0.792     2.208     0.300     0.008 
## trt 
##                1    0.477     0         1           (reference)
##                2    0.523     0.246     1.279     0.203     0.224 
## 
## log(scale)                    2.864    17.537     0.671     0.000 
## log(shape)                    0.075     1.077     0.066     0.261 
## 
## Events                    128 
## Total time at risk         16663 
## Max. log. likelihood      -715.55 
## LR test statistic         65.1 
## Degrees of freedom        8 
## Overall p-value           4.65393e-11
```

Finally, we can also fit a Cox proportional hazard:


```r
model3 <- coxph(y ~ karno + diagtime + age + prior + celltype + trt, 
                data = veteran)
summary(model3)
```

```
## Call:
## coxph(formula = y ~ karno + diagtime + age + prior + celltype + 
##     trt, data = veteran)
## 
##   n= 137, number of events= 128 
## 
##                         coef  exp(coef)   se(coef)      z Pr(>|z|)    
## karno             -3.282e-02  9.677e-01  5.508e-03 -5.958 2.55e-09 ***
## diagtime           8.132e-05  1.000e+00  9.136e-03  0.009  0.99290    
## age               -8.706e-03  9.913e-01  9.300e-03 -0.936  0.34920    
## prior10            7.159e-02  1.074e+00  2.323e-01  0.308  0.75794    
## celltypesquamous  -4.013e-01  6.695e-01  2.827e-01 -1.420  0.15574    
## celltypesmallcell  4.603e-01  1.584e+00  2.662e-01  1.729  0.08383 .  
## celltypeadeno      7.948e-01  2.214e+00  3.029e-01  2.624  0.00869 ** 
## trt2               2.946e-01  1.343e+00  2.075e-01  1.419  0.15577    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##                   exp(coef) exp(-coef) lower .95 upper .95
## karno                0.9677     1.0334    0.9573    0.9782
## diagtime             1.0001     0.9999    0.9823    1.0182
## age                  0.9913     1.0087    0.9734    1.0096
## prior10              1.0742     0.9309    0.6813    1.6937
## celltypesquamous     0.6695     1.4938    0.3847    1.1651
## celltypesmallcell    1.5845     0.6311    0.9403    2.6699
## celltypeadeno        2.2139     0.4517    1.2228    4.0084
## trt2                 1.3426     0.7448    0.8939    2.0166
## 
## Concordance= 0.736  (se = 0.03 )
## Rsquare= 0.364   (max possible= 0.999 )
## Likelihood ratio test= 62.1  on 8 df,   p=1.799e-10
## Wald test            = 62.37  on 8 df,   p=1.596e-10
## Score (logrank) test = 66.74  on 8 df,   p=2.186e-11
```

As we can see, all three models are significant, and they give similar information: ```karno``` and ```celltype``` are significant predictors, both treatment is not.

The method available in this package makes use of *case-base sampling*. That is, person-moments are randomly sampled across the entire follow-up time, with some moments corresponding to cases and others to controls. By sampling person-moments instead of individuals, we can then use logistic regression to fit smooth-in-time parametric hazard functions. See the previous section for more details.

First, we will look at the follow-up time by using population-time plots:


```r
nobs <- nrow(y)
ftime <- veteran$time
ord <- order(ftime, decreasing=TRUE)
plot(0, type='n', xlim=c(0, max(ftime)), ylim=c(0, nobs), 
     xlab='Follow-up time', ylab='Population')
segments(rep(0.0, nobs), 1:nobs, ftime[ord], 1:nobs, col='gray25')
cases <- veteran$status == 1
points((ftime[ord])[cases[ord]], (1:nobs)[cases[ord]], pch=20, col='red', cex=0.5)
```

<img src="smoothHazard_files/figure-html/unnamed-chunk-6-1.png" title="" alt="" style="display: block; margin: auto;" />

Population-time plots are a useful way of visualizing the total follow-up experience, where individuals appear on the y-axis, and follow-up time on the x-axis; each individual's follow-up time is represented by a gray line segment. For convenience, we have ordered the patients according to their time-to-event, and each event is represented by a red dot. The censored observations (of which there is only a few) correspond to the grey lines which do not end with a red dot.

Next, we use case-base sampling to fit a parametric hazard function via logistic regression. First, we will include time as a linear term; as noted above, this corresponds to an Gompertz hazard.


```r
library(casebase)
model4 <- fitSmoothHazard(status ~ time + karno + diagtime + age + prior +
             celltype + trt, data = veteran, ratio=100, type = "uniform")
summary(model4)
```

```
## 
## Call:
## glm(formula = formula, family = binomial(link = link), data = sampleData)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -0.4334  -0.1494  -0.1206  -0.0994   3.4111  
## 
## Coefficients:
##                     Estimate Std. Error z value Pr(>|z|)    
## (Intercept)       -2.7865611  0.7262061  -3.837 0.000124 ***
## time               0.0003885  0.0006458   0.602 0.547478    
## karno             -0.0320120  0.0052834  -6.059 1.37e-09 ***
## diagtime          -0.0007069  0.0093722  -0.075 0.939876    
## age               -0.0050225  0.0092409  -0.544 0.586779    
## prior10            0.0358793  0.2311075   0.155 0.876625    
## celltypesquamous  -0.4118927  0.2835803  -1.452 0.146370    
## celltypesmallcell  0.4395790  0.2624502   1.675 0.093953 .  
## celltypeadeno      0.7470542  0.3007795   2.484 0.013002 *  
## trt2               0.1688823  0.2010431   0.840 0.400891    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 1436.2  on 12927  degrees of freedom
## Residual deviance: 1364.6  on 12918  degrees of freedom
## AIC: 1384.6
## 
## Number of Fisher Scoring iterations: 8
```

Since the output object from ```fitSmoothHazard``` inherits from the ```glm``` class, we see a familiar result when using the function ```summary```.

The main purpose of fitting smooth hazard functions is that it is then relatively easy to compute absolute risks. For example, we can use the function ```absoluteRisk``` to compute the mean absolute risk at 90 days, which can then be compared to the empirical measure.


```r
absoluteRisk(object = model4, time = 90)
```

```
## [1] 0.5778068
```

```r
mean(ftime <= 90)
```

```
## [1] 0.5547445
```

We can also fit a Weibull hazard by using a logarithmic term for time:


```r
model5 <- fitSmoothHazard(status ~ log(time) + karno + diagtime + age + prior +
             celltype + trt, data = veteran, ratio=100, type = "uniform")
summary(model5)
```

```
## 
## Call:
## glm(formula = formula, family = binomial(link = link), data = sampleData)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -0.4267  -0.1497  -0.1190  -0.0974   3.3902  
## 
## Coefficients:
##                    Estimate Std. Error z value Pr(>|z|)    
## (Intercept)       -3.132306   0.749307  -4.180 2.91e-05 ***
## log(time)          0.083588   0.072277   1.156   0.2475    
## karno             -0.031532   0.005447  -5.789 7.10e-09 ***
## diagtime           0.004232   0.009319   0.454   0.6497    
## age               -0.006194   0.009266  -0.669   0.5038    
## prior10           -0.036590   0.229620  -0.159   0.8734    
## celltypesquamous  -0.421632   0.279068  -1.511   0.1308    
## celltypesmallcell  0.495420   0.263323   1.881   0.0599 .  
## celltypeadeno      0.740552   0.301669   2.455   0.0141 *  
## trt2               0.261555   0.203048   1.288   0.1977    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 1436.2  on 12927  degrees of freedom
## Residual deviance: 1366.9  on 12918  degrees of freedom
## AIC: 1386.9
## 
## Number of Fisher Scoring iterations: 8
```

With case-base sampling, it is straightforward to fit a semi-parametric hazard function using splines, which can then be used to estimate the mean absolute risk.


```r
# Fit a spline for time
library(splines)
model6 <- fitSmoothHazard(status ~ bs(time) + karno + diagtime + age + prior +
             celltype + trt, data = veteran, ratio=100, type = "uniform")
summary(model6)
```

```
## 
## Call:
## glm(formula = formula, family = binomial(link = link), data = sampleData)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -0.4523  -0.1540  -0.1194  -0.0962   3.5368  
## 
## Coefficients:
##                    Estimate Std. Error z value Pr(>|z|)    
## (Intercept)       -2.923789   0.718830  -4.067 4.75e-05 ***
## bs(time)1          1.654492   1.022068   1.619  0.10550    
## bs(time)2         -2.786087   1.780166  -1.565  0.11757    
## bs(time)3          1.782795   1.009755   1.766  0.07747 .  
## karno             -0.032825   0.005467  -6.005 1.92e-09 ***
## diagtime          -0.001902   0.009045  -0.210  0.83346    
## age               -0.006261   0.009449  -0.663  0.50760    
## prior10            0.063173   0.234174   0.270  0.78734    
## celltypesquamous  -0.306083   0.281703  -1.087  0.27724    
## celltypesmallcell  0.479477   0.268042   1.789  0.07364 .  
## celltypeadeno      0.874292   0.301055   2.904  0.00368 ** 
## trt2               0.204717   0.204870   0.999  0.31767    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 1436.2  on 12927  degrees of freedom
## Residual deviance: 1365.3  on 12916  degrees of freedom
## AIC: 1389.3
## 
## Number of Fisher Scoring iterations: 8
```

```r
absoluteRisk(object = model6, time = 90)
```

```
## [1] 0.5733611
```

As we can see from the summary, there is little evidence that splines actually improve the fit. Moreover, we can see that estimated individual absolute risks are essentially the same when using either a linear term or splines:


```r
linearRisk <- absoluteRisk(object = model4, time = 90, newdata = veteran)
splineRisk <- absoluteRisk(object = model6, time = 90, newdata = veteran)

plot(linearRisk, splineRisk,
     xlab="Linear", ylab = "Splines", pch=19)
abline(a=0, b=1, lty=2, lwd=2, col='red')
```

<img src="smoothHazard_files/figure-html/unnamed-chunk-11-1.png" title="" alt="" style="display: block; margin: auto;" />

These last three models give similar information as the first three, i.e. the main predictors for the hazard are ```karno``` and ```celltype```, with treatment being non-significant. Moreover, by explicitely including the time variable in the formula, we see that it is not significant; this is evidence that the true hazard is exponential.

<!-- ## Vaccination example -->

<!-- ```{r eval=FALSE} -->
<!-- inpath <- '~/git_repositories/casebase/data/' -->
<!-- etimes <- read.csv(file.path(inpath, 'etimes.csv'), header=TRUE) -->
<!-- vtimes <- read.csv(file.path(inpath, 'vtimes.csv'), header=TRUE) -->
<!-- vtc <- vtimes[etimes[,1],2] -->
<!-- exposed <- (vtc < etimes[,2]) & (vtc >= etimes[,2] - 7.0) -->

<!-- popsize <- nrow(vtimes) -->
<!-- ncases <- nrow(etimes) -->
<!-- n.cases <- rep(1, ncases) -->
<!-- n.exposed.cases <- as.numeric(exposed) -->
<!-- n.unexposed.cases <- as.numeric(!exposed) -->

<!-- # Time-matched sampling: -->

<!-- sample.ctls = function(x) sample(1:popsize, x, FALSE) -->

<!-- ctl.case.ratio = 1 -->
<!-- n.ctls = n.cases * ctl.case.ratio -->
<!-- denom.sample = lapply(n.ctls, sample.ctls) -->
<!-- nctls = sum(n.ctls) -->
<!-- n.exposed.ctls <- rep(NA, ncases) -->

<!-- set.seed(1) -->
<!-- cc <- 0 -->
<!-- casestatus <- NULL -->
<!-- exposure <- NULL -->
<!-- matched <- NULL -->
<!-- for (i in 1:ncases) { -->
<!--     cc <- cc + 1 -->
<!--     exstatus <- (vtimes[denom.sample[[i]],2] < etimes[i,2]) & (vtimes[denom.sample[[i]],2] >= etimes[i,2] - 7.0) -->
<!--     n.exposed.ctls[i] = sum(exstatus) -->
<!--     casestatus <- c(casestatus, rep(1, n.cases[i]), rep(0, n.ctls[i])) -->
<!--     exposure <- c(exposure, rep(0, n.unexposed.cases[i]), rep(1, n.exposed.cases[i]), as.numeric(exstatus)) -->
<!--     matched <- c(matched, rep(cc, n.cases[i] + n.cases[i] * ctl.case.ratio)) -->
<!-- } -->
<!-- # cbind(casestatus, exposure, matched) -->

<!-- n.unexposed.ctls <- n.ctls - n.exposed.ctls -->

<!-- model1 <- clogit(casestatus ~ exposure + strata(as.factor(matched))) -->
<!-- summary(model1) -->

<!-- maxage <- 140/7 -->

<!-- # pdf(file.path(inpath, 'studybase.pdf'), width=8, height=6, paper='special') -->
<!-- postscript(file.path(inpath, 'studybase.eps'), width=8, height=6, paper='special', horizontal=FALSE) -->
<!-- op <- par(mar = c(4.25,4.25,1,1)) -->
<!-- plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population') -->
<!-- rect(0,0,maxage,popsize, col="gray80",border=NA) -->
<!-- legend(0.025 * maxage, 0.975 * popsize, legend=c('Pop.-time'), -->
<!--        pch=c(15), col=c('gray80'), pt.bg=c('gray80'), -->
<!--        cex=1, bg='white', box.col='white') -->
<!-- par(op) -->
<!-- dev.off() -->
<!-- ``` -->

<!-- ```{r eval=FALSE} -->
<!-- idx <- order(vtimes[,2]) -->

<!-- # pdf(file.path(inpath, 'studybaseexposure.pdf'), width=8, height=6, paper='special') -->
<!-- postscript(file.path(inpath, 'studybaseexposure.eps'), width=8, height=6, paper='special', horizontal=FALSE) -->
<!-- op <- par(mar = c(4.25,4.25,1,1)) -->
<!-- plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population') -->
<!-- rect(0,0,maxage,popsize, col="lightblue",border=NA) -->
<!-- segments(vtimes[idx,2]/7, (1:popsize), (vtimes[idx,2]+7)/7, (1:popsize), col='yellow2') -->
<!-- legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.-time','Exposed pop.-time'), -->
<!--        pch=c(15,15), col=c('lightblue','yellow2'), pt.bg=c('lightblue','yellow2'), -->
<!--        cex=1, bg='white', box.col='white') -->
<!-- par(op) -->
<!-- dev.off() -->

<!-- vidx <- order(vtc) -->
<!-- caseidx <- (1:popsize)[vtimes[idx,1] %in% etimes[,1]] -->
<!-- ``` -->

<!-- ```{r eval=FALSE} -->
<!-- # pdf(file.path(inpath, 'caseseries.pdf'), width=8, height=6, paper='special') -->
<!-- postscript(file.path(inpath, 'caseseries.eps'), width=8, height=6, paper='special', horizontal=FALSE) -->
<!-- op <- par(mar = c(4.25,4.25,1,1)) -->
<!-- plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population') -->
<!-- rect(0,0,maxage,popsize, col="lightblue",border=NA) -->
<!-- segments(vtimes[idx,2]/7, (1:popsize), (vtimes[idx,2]+7)/7, (1:popsize), col='yellow2') -->
<!-- points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75) -->
<!-- legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.-time','Exposed pop.-time','Case series'), -->
<!--        pch=c(15,15,20), col=c('lightblue','yellow2','red'), pt.bg=c('lightblue','yellow2','red'), -->
<!--        cex=1, bg='white', box.col='white') -->
<!-- par(op) -->
<!-- dev.off() -->
<!-- ``` -->

<!-- ```{r eval=FALSE} -->
<!-- # pdf(file.path(inpath, 'baseseries.pdf'), width=8, height=6, paper='special') -->
<!-- postscript(file.path(inpath, 'baseseries.eps'), width=8, height=6, paper='special', horizontal=FALSE) -->
<!-- op <- par(mar = c(4.25,4.25,1,1)) -->
<!-- plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population') -->
<!-- rect(0,0,maxage,popsize, col="lightblue",border=NA) -->
<!-- segments(vtimes[idx,2]/7, (1:popsize), (vtimes[idx,2]+7)/7, (1:popsize), col='yellow2') -->
<!-- points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75) -->
<!-- points(etimes[matched[casestatus==0],2]/7, vtimes[unlist(denom.sample),1], pch=20, col='blue', cex=0.75) -->
<!-- legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.-time','Exposed pop.-time','Case series','Base series'), -->
<!--        pch=c(15,15,20,20), col=c('lightblue','yellow2','red','blue'), pt.bg=c('lightblue','yellow2','red','blue'), -->
<!--        cex=1, bg='white', box.col='white') -->
<!-- par(op) -->
<!-- dev.off() -->
<!-- ``` -->

<!-- ```{r eval=FALSE} -->
<!-- # Self-matched sampling: -->

<!-- minday <- 0 -->
<!-- maxday <- 140 -->

<!-- nbase <- 1 -->
<!-- times <- rep(NA, 2 * (ncases + nbase * ncases)) -->
<!-- casestatus <- rep(NA, 2 * (ncases + nbase * ncases)) -->
<!-- exposure <- rep(NA, 2 * (ncases + nbase * ncases)) -->
<!-- matched <- rep(NA, 2 * (ncases + nbase * ncases)) -->

<!-- set.seed(1) -->
<!-- counter <- 0 -->
<!-- casecounter <- 0 -->
<!-- for (i in 1:ncases) { -->
<!--     if (is.finite(etimes[i,2]) & etimes[i,2] > minday & etimes[i,2] <= maxday) { -->
<!--         counter <- counter + 1 -->
<!--         casecounter <- casecounter + 1 -->
<!--         times[counter] <- etimes[i,2] -->
<!--         casestatus[counter] <- 1 -->
<!--         exposure[counter] <- (etimes[i,2] - vtc[i]) > 0.0 & (etimes[i,2] - vtc[i]) < 7.0 -->
<!--         matched[counter] <- casecounter -->
<!--         # n <- rpois(1, nbase) -->
<!--         n <- nbase -->
<!--         for (j in 1:n) { -->
<!--             counter <- counter + 1 -->
<!--             u <- runif(1, 0.0, nbase) -->
<!--             times[counter] <- runif(1, min=minday, max=maxday) -->
<!--             casestatus[counter] <- 0 -->
<!--             exposure[counter] <- (times[counter] - vtc[i]) > 0.0 & (times[counter] - vtc[i]) < 7.0 -->
<!--             matched[counter] <- casecounter -->
<!--         } -->
<!--     } -->
<!--     # if (i %% 1000 == 0) -->
<!--     #     print(i) -->
<!-- } -->
<!-- # counter -->
<!-- # casecounter -->
<!-- times <- times[1:counter] -->
<!-- casestatus <- casestatus[1:counter] -->
<!-- exposure <- exposure[1:counter] -->
<!-- matched <- matched[1:counter] -->

<!-- times2 <- times^2 -->
<!-- model2 <- clogit(casestatus ~ exposure + times + times2 + strata(as.factor(matched))) -->
<!-- summary(model2) -->

<!-- model3 <- clogit(casestatus ~ exposure + bs(times) + strata(as.factor(matched))) -->
<!-- summary(model3) -->
<!-- ``` -->

<!-- ```{r eval=FALSE} -->
<!-- # pdf(file.path(inpath, 'caseseries2.pdf'), width=8, height=6, paper='special') -->
<!-- postscript(file.path(inpath, 'caseseries2.eps'), width=8, height=6, paper='special', horizontal=FALSE) -->
<!-- op <- par(mar = c(4.25,4.25,1,1)) -->
<!-- plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population') -->
<!-- rect(0,0,maxage,popsize, col="lightblue",border=NA) -->
<!-- segments(vtimes[idx,2]/7, (1:popsize), (vtimes[idx,2]+7)/7, (1:popsize), col='yellow2') -->
<!-- points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75) -->
<!-- legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.n-time','Exposed pop.-time','Case series'), -->
<!--        pch=c(15,15,20), col=c('lightblue','yellow2','red'), pt.bg=c('lightblue','yellow2','red'), -->
<!--        cex=1, bg='white', box.col='white') -->
<!-- par(op) -->
<!-- dev.off() -->
<!-- ``` -->

<!-- ```{r eval=FALSE} -->
<!-- # pdf(file.path(inpath, 'selfmatched.pdf'), width=8, height=6, paper='special') -->
<!-- postscript(file.path(inpath, 'selfmatched.eps'), width=8, height=6, paper='special', horizontal=FALSE) -->
<!-- op <- par(mar = c(4.25,4.25,1,1)) -->
<!-- plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population') -->
<!-- rect(0,0,maxage,popsize, col="white",border=NA) -->
<!-- segments(rep(minday/7, ncases), caseidx, rep(maxday/7, ncases), caseidx, col='lightblue') -->
<!-- segments(vtc[vidx]/7, caseidx, (vtc[vidx]+7)/7, caseidx, col='yellow2') -->
<!-- points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75) -->
<!-- legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.-time','Exposed pop.-time','Case series'), -->
<!--        pch=c(15,15,20), col=c('lightblue','yellow2','red'), pt.bg=c('lightblue','yellow2','red'), -->
<!--        cex=1, bg='white', box.col='white') -->
<!-- par(op) -->
<!-- dev.off() -->
<!-- ``` -->

<!-- ```{r eval=FALSE} -->
<!-- # pdf(file.path(inpath, 'baseseries2.pdf'), width=8, height=6, paper='special') -->
<!-- postscript(file.path(inpath, 'baseseries2.eps'), width=8, height=6, paper='special', horizontal=FALSE) -->
<!-- op <- par(mar = c(4.25,4.25,1,1)) -->
<!-- plot(1, 1, xlim=c(0,maxage), ylim=c(0,popsize), type='n', xlab='Age (weeks)', ylab='Population') -->
<!-- rect(0,0,maxage,popsize, col="white",border=NA) -->
<!-- segments(rep(minday/7, ncases), caseidx, rep(maxday/7, ncases), caseidx, col='lightblue') -->
<!-- segments(vtc[vidx]/7, caseidx, (vtc[vidx]+7)/7, caseidx, col='yellow2') -->
<!-- points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75) -->
<!-- points(times[casestatus==0]/7, caseidx[matched[casestatus==0]], pch=20, col='blue', cex=0.75) -->
<!-- legend(0.025 * maxage, 0.975 * popsize, legend=c('Unexposed pop.-time','Exposed pop.-time','Case series','Base series'), -->
<!--        pch=c(15,15,20,20), col=c('lightblue','yellow2','red','blue'), pt.bg=c('lightblue','yellow2','red','blue'), -->
<!--        cex=1, bg='white', box.col='white') -->
<!-- par(op) -->
<!-- dev.off() -->
<!-- ``` -->

<!-- ```{r eval=FALSE} -->
<!-- # pdf(file.path(inpath, 'baseseries3.pdf'), width=8, height=6, paper='special') -->
<!-- postscript(file.path(inpath, 'baseseries3.eps'), width=8, height=6, paper='special', horizontal=FALSE) -->
<!-- op <- par(mar = c(4.25,4.25,1,1)) -->
<!-- plot(1, 1, xlim=c(0,maxage), ylim=c(0,1000), type='n', xlab='Age (weeks)', ylab='Population') -->
<!-- rect(0,0,maxage,popsize, col="white",border=NA) -->
<!-- segments(rep(minday/7, ncases), caseidx, rep(maxday/7, ncases), caseidx, col='lightblue') -->
<!-- segments(vtc[vidx]/7, caseidx, (vtc[vidx]+7)/7, caseidx, col='yellow2') -->
<!-- points(etimes[vidx,2]/7, caseidx, pch=20, col='red', cex=0.75) -->
<!-- points(times[casestatus==0]/7, caseidx[matched[casestatus==0]], pch=20, col='blue', cex=0.75) -->
<!-- par(op) -->
<!-- dev.off() -->
<!-- ``` -->

## Session information


```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.4 LTS
## 
## locale:
##  [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8    
##  [5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
##  [7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] splines   stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
## [1] eha_2.4-3         survival_2.38-3   VGAM_1.0-0        casebase_0.0.9000
## 
## loaded via a namespace (and not attached):
##  [1] digest_0.6.9    formatR_1.2.1   magrittr_1.5    evaluate_0.8   
##  [5] stringi_1.0-1   rmarkdown_0.9.2 devtools_1.10.0 tools_3.2.3    
##  [9] stringr_1.0.0   yaml_2.1.13     memoise_1.0.0   htmltools_0.3  
## [13] knitr_1.12.3
```
## References
