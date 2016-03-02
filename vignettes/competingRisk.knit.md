---
title: "Competing risk analysis using case-base sampling"
author: "Maxime Turgeon"
date: "2016-02-29"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: reference.bib
---

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

## Population-time plot

In order to try and visualize the incidence density of each event, we can look at a population-time plot: on the Y axis, we order the failure times from shortest (at the top) to longest (at the bottom). Then each line corresponds to the follow-up time of one individual. Failure times associated to the event of interest or a competing event can then be highlighted on the plot using coloured dots.


```r
nobs <- nrow(DT)
ftime <- DT$ftime
ord <- order(ftime, decreasing=TRUE)
plot(0, type='n', xlim=c(0, max(ftime)), ylim=c(0, nobs), 
     xlab='Follow-up time', ylab='Population')
#segments(rep(0.0, nobs), 1:nobs, ftime[ord], 1:nobs, col='gray25')
polygon(c(0, max(ftime), ftime[ord], 0), c(0, 0, 1:nobs, nobs), col = "gray90")
cases <- DT$Status %in% c(1, 2)
colour <- c("red", "blue")[DT$Status[cases]]

# randomly move the cases vertically
moved_cases <- sapply((1:nobs)[cases[ord]], sample, size=1)
points((ftime[ord])[cases[ord]], moved_cases, pch=20, col=colour, cex=0.5)
legend("topright", legend=c("Relapse", "Competing event"), col=c("red", "blue"),
       pch=20)
```

<img src="competingRisk_files/figure-html/unnamed-chunk-2-1.png" title="" alt="" style="display: block; margin: auto;" />

## Analysis


```r
library(casebase)
library(VGAM)
model1 <- fitSmoothHazard(Status ~ ftime + Sex + D + Phase + Source + Age, 
                          data = DT, ratio=100, type = "uniform", time="ftime")
summary(model1)
```

```
## 
## Call:
## VGAM::vglm(formula = formula, family = VGAM::multinomial, data = combData)
## 
## Pearson residuals:
##                        Min       1Q   Median        3Q   Max
## log(mu[,1]/mu[,3]) -0.2261 -0.07052 -0.03950 -0.014747 46.46
## log(mu[,2]/mu[,3]) -0.3266 -0.09127 -0.03727 -0.008339 20.05
## 
## Coefficients:
##                Estimate Std. Error z value Pr(>|z|)    
## (Intercept):1  -3.36554    0.68395  -4.921 8.62e-07 ***
## (Intercept):2  -2.52403    0.46645  -5.411 6.26e-08 ***
## ftime:1        -0.06967    0.01478  -4.714 2.43e-06 ***
## ftime:2        -0.10338    0.01825  -5.665 1.47e-08 ***
## SexM:1         -0.25845    0.28272  -0.914 0.360625    
## SexM:2         -0.35062    0.23646  -1.483 0.138124    
## DAML:1         -0.58347    0.30123  -1.937 0.052755 .  
## DAML:2         -0.08431    0.27826  -0.303 0.761888    
## PhaseCR2:1      0.08166    0.46576   0.175 0.860821    
## PhaseCR2:2      0.18169    0.33029   0.550 0.582270    
## PhaseCR3:1      0.36205    0.69075   0.524 0.600180    
## PhaseCR3:2      0.10690    0.52601   0.203 0.838952    
## PhaseRelapse:1  1.37533    0.39212   3.507 0.000452 ***
## PhaseRelapse:2  0.70390    0.30855   2.281 0.022531 *  
## SourcePB:1      0.38651    0.57251   0.675 0.499606    
## SourcePB:2     -1.04172    0.35633  -2.923 0.003461 ** 
## Age:1          -0.00758    0.01196  -0.634 0.526026    
## Age:2           0.02588    0.01003   2.580 0.009883 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Number of linear predictors:  2 
## 
## Names of linear predictors: log(mu[,1]/mu[,3]), log(mu[,2]/mu[,3])
## 
## Dispersion Parameter for multinomial family:   1
## 
## Residual deviance: 1416.047 on 26444 degrees of freedom
## 
## Log-likelihood: -708.0236 on 26444 degrees of freedom
## 
## Number of iterations: 10
```


```r
model2 <- fitSmoothHazard(Status ~ log(ftime) + Sex + D + Phase + Source + Age, 
                          data = DT, ratio=100, type = "uniform", time="ftime")
summary(model2)
```

```
## 
## Call:
## VGAM::vglm(formula = formula, family = VGAM::multinomial, data = combData)
## 
## Pearson residuals:
##                        Min       1Q   Median       3Q   Max
## log(mu[,1]/mu[,3]) -0.4364 -0.06846 -0.04617 -0.03488 30.05
## log(mu[,2]/mu[,3]) -0.4798 -0.07737 -0.05445 -0.04416 21.69
## 
## Coefficients:
##                 Estimate Std. Error z value Pr(>|z|)    
## (Intercept):1  -3.886137   0.703454  -5.524 3.31e-08 ***
## (Intercept):2  -3.040376   0.462164  -6.579 4.75e-11 ***
## log(ftime):1   -0.348206   0.070358  -4.949 7.46e-07 ***
## log(ftime):2   -0.428307   0.057147  -7.495 6.64e-14 ***
## SexM:1         -0.445740   0.295271  -1.510 0.131146    
## SexM:2         -0.503866   0.241919  -2.083 0.037271 *  
## DAML:1         -0.716902   0.304968  -2.351 0.018736 *  
## DAML:2         -0.231844   0.286188  -0.810 0.417876    
## PhaseCR2:1      0.234525   0.470276   0.499 0.617993    
## PhaseCR2:2      0.356388   0.333476   1.069 0.285201    
## PhaseCR3:1      0.464975   0.708295   0.656 0.511521    
## PhaseCR3:2      0.138775   0.531549   0.261 0.794034    
## PhaseRelapse:1  1.503120   0.394482   3.810 0.000139 ***
## PhaseRelapse:2  0.893981   0.310009   2.884 0.003930 ** 
## SourcePB:1      0.633374   0.600817   1.054 0.291797    
## SourcePB:2     -0.968908   0.367066  -2.640 0.008300 ** 
## Age:1          -0.003343   0.011817  -0.283 0.777242    
## Age:2           0.030514   0.010055   3.035 0.002406 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Number of linear predictors:  2 
## 
## Names of linear predictors: log(mu[,1]/mu[,3]), log(mu[,2]/mu[,3])
## 
## Dispersion Parameter for multinomial family:   1
## 
## Residual deviance: 1501.208 on 26444 degrees of freedom
## 
## Log-likelihood: -750.604 on 26444 degrees of freedom
## 
## Number of iterations: 8
```


```r
linearRisk <- absoluteRisk(object = model1, time = 60, method = "montecarlo")
logRisk <- absoluteRisk(object = model2, time = 60, method = "montecarlo")

plot(linearRisk[,1], logRisk[,1],
     xlab="Linear", ylab = "Log", pch=19, cex=0.5)
points(linearRisk[,2], logRisk[,2],
       col = 'blue', pch=19, cex=0.5)
abline(a=0, b=1, lty=2, lwd=2, col='red')
```

<img src="competingRisk_files/figure-html/unnamed-chunk-5-1.png" title="" alt="" style="display: block; margin: auto;" />

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
## [1] VGAM_1.0-0        casebase_0.0.9000
## 
## loaded via a namespace (and not attached):
##  [1] digest_0.6.9    formatR_1.2.1   magrittr_1.5    evaluate_0.8   
##  [5] stringi_1.0-1   rmarkdown_0.9.2 devtools_1.10.0 tools_3.2.3    
##  [9] stringr_1.0.0   yaml_2.1.13     survival_2.38-3 memoise_1.0.0  
## [13] htmltools_0.3   knitr_1.12.3
```
## References
