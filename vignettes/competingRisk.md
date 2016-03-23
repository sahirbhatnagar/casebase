# Competing risk analysis using case-base sampling
Maxime Turgeon  
`r Sys.Date()`  

## Data

We will use the same data that was used in Scrucca *et al* [-@scrucca2010regression]. The data is available on the main author's [website](http://www.stat.unipg.it/luca/R/).


```r
DT <- read.csv(system.file("extdata", "bmtcrr.csv", package = "casebase"))
head(DT)
```

```
##   Sex   D   Phase Age Status Source  ftime
## 1   M ALL Relapse  48      2  BM+PB   0.67
## 2   F AML     CR2  23      1  BM+PB   9.50
## 3   M ALL     CR3   7      0  BM+PB 131.77
## 4   F ALL     CR2  26      2  BM+PB  24.03
## 5   F ALL     CR2  36      2  BM+PB   1.47
## 6   M ALL Relapse  17      2  BM+PB   2.23
```

We will perform a competing risk analysis on data from 177 patients who received a stem cell transplant for acute leukemia. The event of interest in relapse, but other competing causes (e.g. transplant-related death) need to be taken into account. We also want to take into account the effect of several covariates such as Sex, Disease (lymphoblastic or myeloblastic leukemia, abbreviated as ALL and AML, respectively), Phase at transplant (Relapse, CR1, CR2, CR3), Source of stem cells (bone marrow and peripheral blood, coded as BM+PB, or peripheral blood, coded as PB), and Age. Below, we reproduce their Table 1:

| Variable | Description            | Statistical summary    |
| -------- | ---------------------- | ---------------------- |
| Sex      | Sex                    | M=Male (100) <br> F=Female (77) |
| D        | Disease                | ALL (73) <br> AML (104) |
| Phase    | Phase                  | CR1 (47) <br> CR2 (45) <br> CR3 (12) <br> Relapse (73) |
| Source   | Type of transplant     | BM+PB (21) <br> PB (156) |
| Age      | Age of patient (years) | 4–62 <br> 30.47 (13.04) |
| Ftime    | Failure time (months)  | 0.13–131.77 <br> 20.28 (30.78) |
| Status   | Status indicator       | 0=censored (46) <br> 1=relapse (56) <br> 2=competing event (75) |

The statistical summary is generated differently for continuous and categorical variables: 

 - For continuous variables, we are given the range, followed by the mean and standard deviation.
 
 - For categorical variables, we are given the counts for each category.
 
Note that failure time can also correspond to censoring.

## Population-time plots

In order to try and visualize the incidence density of relapse, we can look at a population-time plot: on the X-axis we have time, and on the Y-axis we have the size of the risk set at a particular time point. Failure times associated to the event of interest can then be highlighted on the plot using red dots.


```r
nobs <- nrow(DT)
ftime <- DT$ftime
ord <- order(ftime, decreasing=FALSE)

# We split the person-moments in four categories:
# 1) at-risk
# 2) main event
# 3) competing event
# 4) censored
yCoords <- cbind(cumsum(DT[ord, "Status"] == 2), 
                 cumsum(DT[ord, "Status"] == 1),
                 cumsum(DT[ord, "Status"] == 0))
yCoords <- cbind(yCoords, nobs - rowSums(yCoords))

# Plot only at-risk
plot(0, type='n', xlim=c(0, max(ftime)), ylim=c(0, nobs), 
     xlab='Follow-up time', ylab='Population')
polygon(c(0, 0, ftime[ord], max(ftime), 0),
        c(0, nobs, yCoords[,4], 0, 0), col = "grey90")
cases <- DT[, "Status"] == 1

# randomly move the cases vertically
moved_cases <- yCoords[cases[ord], 4] * runif(sum(cases))
points((ftime[ord])[cases[ord]], moved_cases, pch=20, col="red", cex=1)
```

<img src="competingRisk_files/figure-html/poptime1-1.png" title="" alt="" style="display: block; margin: auto;" />

We can right away draw a few conclusions from this plot: first of all, we get a sense of how quickly the size of the risk set changes over time. We also see that the incidence density is non-constant: most relapses occur before 15 months. Finally, we also see that the risk set keeps shrinking after the last event has occured; this could be due to either censoring or the competing event.

To get an idea of whether only relapse is responsible for the shrinking of the risk set in the first few months of follow-up, we can also keep track of how many events have occured at each time point:


```r
# Plot at-risk and events
plot(0, type='n', xlim=c(0, max(ftime)), ylim=c(0, nobs), 
     xlab='Follow-up time', ylab='Population')
polygon(c(0,ftime[ord], max(ftime), 0), c(0, yCoords[,2], 0, 0), col = "firebrick3")
polygon(c(0, ftime[ord], ftime[rev(ord)], 0, 0),
        c(0, yCoords[,2], rev(yCoords[,2] + yCoords[,4]), nobs, 0), col = "grey90")

# randomly move the cases vertically
moved_cases <- yCoords[cases[ord], 2] + yCoords[cases[ord], 4] * runif(sum(cases))
points((ftime[ord])[cases[ord]], moved_cases, pch=20, col="red", cex=1)
legend("topright", legend=c("Relapse", "At-risk"), 
       col=c("firebrick3", "grey90"),
       pch=15)
```

<img src="competingRisk_files/figure-html/poptime2-1.png" title="" alt="" style="display: block; margin: auto;" />

Therefore, there is also censoring and loss due to competing events happening in the first few months. However, with this plot, we can't differentiate bwetween the two contributions. For this reason we can also keep track of the number of competing events at each time point:


```r
plot(0, type='n', xlim=c(0, max(ftime)), ylim=c(0, nobs), 
     xlab='Follow-up time', ylab='Population')
polygon(c(0, max(ftime), max(ftime), 0),
        c(0, 0, nobs, nobs), col = "white")
# Event of interest
polygon(c(0,ftime[ord], max(ftime), 0), c(0, yCoords[,2], 0, 0), col = "firebrick3")
# Risk set
polygon(c(0, ftime[ord], ftime[rev(ord)], 0, 0),
        c(0, yCoords[,2], rev(yCoords[,2] + yCoords[,4]), nobs, 0), col = "grey90")
# Competing event
polygon(c(0, ftime[ord], max(ftime), 0), c(nobs, nobs - yCoords[,1], nobs, nobs), col = "dodgerblue2")

# randomly move the cases vertically
moved_cases <- yCoords[cases[ord], 2] + yCoords[cases[ord], 4] * runif(sum(cases))
points((ftime[ord])[cases[ord]], moved_cases, pch=20, col="red", cex=1)
legend("topright", legend=c("Relapse", "Competing event", "At-risk"), 
       col=c("firebrick3", "dodgerblue2", "grey90"),
       pch=15)
```

<img src="competingRisk_files/figure-html/poptime3-1.png" title="" alt="" style="display: block; margin: auto;" />

From this last plot, we can see that there is no censoring during the first 10 months. Moreover, we see that the last competing event occurs around 20 months. Putting all this information together, we have evidence of two types of patients: very sick patients who either relapse or have a competing event early on, and healthier patients who are eventually lost to follow-up.

## Analysis

We now turn to the analysis of this dataset. The population-time plots above give evidence of non-constant hazard; therefore, we will explicitely include time in the model. Note that we also include all other variables as possible confounders. First, we include time as a linear term:


```r
library(casebase)
model1 <- fitSmoothHazard(Status ~ ftime + Sex + D + Phase + Source + Age, 
                          data = DT, ratio=1000, type = "uniform", time="ftime")
summary(model1)
```

```
## 
## Call:
## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
## 
## Pearson residuals:
##                         Min       1Q   Median        3Q    Max
## log(mu[,2]/mu[,1]) -0.07345 -0.02205 -0.01239 -0.004704 144.73
## log(mu[,3]/mu[,1]) -0.10895 -0.02841 -0.01179 -0.002704  63.41
## 
## Coefficients:
##                 Estimate Std. Error z value Pr(>|z|)    
## (Intercept):1  -3.474588   0.680663  -5.105 3.31e-07 ***
## (Intercept):2  -2.549614   0.459171  -5.553 2.81e-08 ***
## ftime:1        -0.069401   0.014631  -4.744 2.10e-06 ***
## ftime:2        -0.103346   0.018076  -5.717 1.08e-08 ***
## SexM:1         -0.293579   0.281057  -1.045  0.29623    
## SexM:2         -0.408156   0.234223  -1.743  0.08140 .  
## DAML:1         -0.629272   0.298733  -2.106  0.03516 *  
## DAML:2         -0.136885   0.274495  -0.499  0.61800    
## PhaseCR2:1      0.191556   0.465535   0.411  0.68072    
## PhaseCR2:2      0.290034   0.329955   0.879  0.37940    
## PhaseCR3:1      0.509231   0.690816   0.737  0.46104    
## PhaseCR3:2      0.248652   0.523925   0.475  0.63508    
## PhaseRelapse:1  1.452848   0.390601   3.720  0.00020 ***
## PhaseRelapse:2  0.779157   0.306387   2.543  0.01099 *  
## SourcePB:1      0.454649   0.568322   0.800  0.42372    
## SourcePB:2     -1.079430   0.353163  -3.056  0.00224 ** 
## Age:1          -0.006372   0.011866  -0.537  0.59125    
## Age:2           0.028030   0.009902   2.831  0.00464 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Number of linear predictors:  2 
## 
## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
## 
## Dispersion Parameter for multinomial family:   1
## 
## Residual deviance: 2009.985 on 262244 degrees of freedom
## 
## Log-likelihood: -1004.992 on 262244 degrees of freedom
## 
## Number of iterations: 13 
## 
## Reference group is level  1  of the response
```

Because of the results in Turgeon *et al* [-@turgeonCompRisk], the standard errors we obtain from the multinomial logit fit are asymptotically correct, and therefore can be used to construct asymptotic confidence intervals. 

From this summary, we see that time is indeed significant, as is Phase (only relapse vs. CR1). Interestingly, we see that the type of disease is only significant for the event of interest, whereas the type of transplant and the age of the patient are only significant for the competing event.

Next, we include the logarithm of time in the model (which leads to a Weibull hazard):



```r
model2 <- fitSmoothHazard(Status ~ log(ftime) + Sex + D + Phase + Source + Age, 
                          data = DT, ratio=1000, type = "uniform", time="ftime")
summary(model2)
```

```
## 
## Call:
## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
## 
## Pearson residuals:
##                        Min       1Q   Median       3Q   Max
## log(mu[,2]/mu[,1]) -0.1723 -0.02169 -0.01484 -0.01120 91.45
## log(mu[,3]/mu[,1]) -0.2699 -0.02450 -0.01754 -0.01444 69.01
## 
## Coefficients:
##                 Estimate Std. Error z value Pr(>|z|)    
## (Intercept):1  -3.931793   0.701633  -5.604 2.10e-08 ***
## (Intercept):2  -3.043980   0.459504  -6.624 3.48e-11 ***
## log(ftime):1   -0.326174   0.067934  -4.801 1.58e-06 ***
## log(ftime):2   -0.405412   0.054786  -7.400 1.36e-13 ***
## SexM:1         -0.452946   0.291820  -1.552 0.120628    
## SexM:2         -0.515177   0.239316  -2.153 0.031342 *  
## DAML:1         -0.682449   0.303771  -2.247 0.024666 *  
## DAML:2         -0.149322   0.287263  -0.520 0.603197    
## PhaseCR2:1      0.261893   0.467302   0.560 0.575181    
## PhaseCR2:2      0.385421   0.329610   1.169 0.242273    
## PhaseCR3:1      0.438951   0.710562   0.618 0.536738    
## PhaseCR3:2      0.082495   0.532621   0.155 0.876912    
## PhaseRelapse:1  1.497704   0.392762   3.813 0.000137 ***
## PhaseRelapse:2  0.880194   0.307423   2.863 0.004195 ** 
## SourcePB:1      0.651808   0.600985   1.085 0.278114    
## SourcePB:2     -0.971796   0.368542  -2.637 0.008368 ** 
## Age:1          -0.004476   0.011733  -0.382 0.702822    
## Age:2           0.027958   0.009871   2.832 0.004620 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Number of linear predictors:  2 
## 
## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
## 
## Dispersion Parameter for multinomial family:   1
## 
## Residual deviance: 2104.931 on 262244 degrees of freedom
## 
## Log-likelihood: -1052.465 on 262244 degrees of freedom
## 
## Number of iterations: 11 
## 
## Reference group is level  1  of the response
```

As we can see, the results are similar to the ones with a Gompertz hazard, although Sex is now significant for the competing event.

Finally, using splines, we can be quite flexible about the way the hazard depends on time:


```r
model3 <- fitSmoothHazard(Status ~ bs(ftime) + Sex + D + Phase + Source + Age, 
                          data = DT, ratio=1000, type = "uniform", time="ftime")
summary(model3)
```

```
## 
## Call:
## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
## 
## Pearson residuals:
##                         Min       1Q    Median         3Q   Max
## log(mu[,2]/mu[,1]) -0.06477 -0.02234 -0.012758 -2.412e-03 177.3
## log(mu[,3]/mu[,1]) -0.08998 -0.03014 -0.004554 -2.306e-06 118.6
## 
## Coefficients:
##                  Estimate Std. Error z value Pr(>|z|)    
## (Intercept):1   -3.741211   0.704544  -5.310 1.10e-07 ***
## (Intercept):2   -3.269145   0.506820  -6.450 1.12e-10 ***
## bs(ftime)1:1     0.007116   2.263182   0.003 0.997491    
## bs(ftime)1:2     7.001805   3.651896   1.917 0.055199 .  
## bs(ftime)2:1   -16.189451   8.163262  -1.983 0.047344 *  
## bs(ftime)2:2   -77.534046  25.577967  -3.031 0.002435 ** 
## bs(ftime)3:1    -2.417463  10.211989  -0.237 0.812868    
## bs(ftime)3:2    -2.437630  22.304549  -0.109 0.912974    
## SexM:1          -0.287843   0.282408  -1.019 0.308086    
## SexM:2          -0.402536   0.235158  -1.712 0.086939 .  
## DAML:1          -0.632579   0.299835  -2.110 0.034879 *  
## DAML:2          -0.163313   0.273080  -0.598 0.549813    
## PhaseCR2:1       0.159591   0.465919   0.343 0.731953    
## PhaseCR2:2       0.292454   0.330122   0.886 0.375673    
## PhaseCR3:1       0.523730   0.693387   0.755 0.450056    
## PhaseCR3:2       0.261503   0.524775   0.498 0.618262    
## PhaseRelapse:1   1.476824   0.394558   3.743 0.000182 ***
## PhaseRelapse:2   0.877603   0.310625   2.825 0.004724 ** 
## SourcePB:1       0.431423   0.572585   0.753 0.451170    
## SourcePB:2      -1.131722   0.356749  -3.172 0.001512 ** 
## Age:1           -0.004786   0.012043  -0.397 0.691044    
## Age:2            0.030841   0.010086   3.058 0.002229 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Number of linear predictors:  2 
## 
## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
## 
## Dispersion Parameter for multinomial family:   1
## 
## Residual deviance: 1995.605 on 262240 degrees of freedom
## 
## Log-likelihood: -997.8027 on 262240 degrees of freedom
## 
## Number of iterations: 18 
## 
## Reference group is level  1  of the response
```

Again, we see that the results are quite similar for this third model.

### Absolute risk

We now look at the 2-year risk of relapse:


```r
linearRisk <- absoluteRisk(object = model1, time = 24, newdata = DT[1:10,])
logRisk <- absoluteRisk(object = model2, time = 24, newdata = DT[1:10,])
splineRisk <- absoluteRisk(object = model3, time = 24, newdata = DT[1:10,])
```


```r
plot(linearRisk[,1], logRisk[,1],
     xlab="Linear", ylab = "Log/Spline", pch=19,
     xlim=c(0,1), ylim=c(0,1), col='red')
points(linearRisk[,1], splineRisk[,1],
       col = 'blue', pch=19)
abline(a=0, b=1, lty=2, lwd=2)
legend("topleft", legend=c("Log", "Spline"),
       pch=19, col=c("red", "blue"))
```

<img src="competingRisk_files/figure-html/absRiskPlot-1.png" title="" alt="" style="display: block; margin: auto;" />

As we can see, Model 1 and Model 2 give different absolute risk predictions, but the linear and the spline model actually give very similar results. We can also estimate the mean absolute risk for the entire dataset:

```r
# The first column corresponds to the event of interest
mean(linearRisk[,1])
```

```
## [1] 0.1958714
```

```r
mean(logRisk[,1])
```

```
## [1] 0.1597606
```

```r
mean(splineRisk[,1])
```

```
## [1] 0.1970028
```

## Session information


```
## R version 3.2.4 Revised (2016-03-16 r70336)
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
## [1] casebase_0.0.9000 VGAM_1.0-1       
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.3      knitr_1.12.3     magrittr_1.5     devtools_1.10.0 
##  [5] munsell_0.4.2    colorspace_1.2-6 R6_2.1.2         httr_1.1.0      
##  [9] stringr_1.0.0    plyr_1.8.3       tools_3.2.4      grid_3.2.4      
## [13] data.table_1.9.6 gtable_0.1.2     withr_1.0.1      git2r_0.13.1    
## [17] htmltools_0.3    yaml_2.1.13      survival_2.38-3  digest_0.6.9    
## [21] ggplot2_2.0.0    formatR_1.3      codetools_0.2-14 curl_0.9.5      
## [25] memoise_1.0.0    evaluate_0.8.3   rmarkdown_0.9.5  stringi_1.0-1   
## [29] scales_0.3.0     chron_2.3-47
```
## References
