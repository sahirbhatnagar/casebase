Data
----

We will use the same data that was used in Scrucca *et al* \[-@scrucca2010regression\]. The data is available on the main author's [website](http://www.stat.unipg.it/luca/R/).

``` r
DT <- read.csv(system.file("extdata", "bmtcrr.csv", package = "casebase"))
head(DT)
```

    ##   Sex   D   Phase Age Status Source  ftime
    ## 1   M ALL Relapse  48      2  BM+PB   0.67
    ## 2   F AML     CR2  23      1  BM+PB   9.50
    ## 3   M ALL     CR3   7      0  BM+PB 131.77
    ## 4   F ALL     CR2  26      2  BM+PB  24.03
    ## 5   F ALL     CR2  36      2  BM+PB   1.47
    ## 6   M ALL Relapse  17      2  BM+PB   2.23

We will perform a competing risk analysis on data from 177 patients who received a stem cell transplant for acute leukemia. The event of interest in relapse, but other competing causes (e.g. transplant-related death) need to be taken into account. We also want to take into account the effect of several covariates such as Sex, Disease (lymphoblastic or myeloblastic leukemia, abbreviated as ALL and AML, respectively), Phase at transplant (Relapse, CR1, CR2, CR3), Source of stem cells (bone marrow and peripheral blood, coded as BM+PB, or peripheral blood, coded as PB), and Age. Below, we reproduce their Table 1:

<table style="width:76%;">
<colgroup>
<col width="12%" />
<col width="31%" />
<col width="31%" />
</colgroup>
<thead>
<tr class="header">
<th>Variable</th>
<th>Description</th>
<th>Statistical summary</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Sex</td>
<td>Sex</td>
<td>M=Male (100) <br> F=Female (77)</td>
</tr>
<tr class="even">
<td>D</td>
<td>Disease</td>
<td>ALL (73) <br> AML (104)</td>
</tr>
<tr class="odd">
<td>Phase</td>
<td>Phase</td>
<td>CR1 (47) <br> CR2 (45) <br> CR3 (12) <br> Relapse (73)</td>
</tr>
<tr class="even">
<td>Source</td>
<td>Type of transplant</td>
<td>BM+PB (21) <br> PB (156)</td>
</tr>
<tr class="odd">
<td>Age</td>
<td>Age of patient (years)</td>
<td>4–62 <br> 30.47 (13.04)</td>
</tr>
<tr class="even">
<td>Ftime</td>
<td>Failure time (months)</td>
<td>0.13–131.77 <br> 20.28 (30.78)</td>
</tr>
<tr class="odd">
<td>Status</td>
<td>Status indicator</td>
<td>0=censored (46) <br> 1=relapse (56) <br> 2=competing event (75)</td>
</tr>
</tbody>
</table>

The statistical summary is generated differently for continuous and categorical variables:

-   For continuous variables, we are given the range, followed by the mean and standard deviation.

-   For categorical variables, we are given the counts for each category.

Note that failure time can also correspond to censoring.

Population-time plots
---------------------

In order to try and visualize the incidence density of relapse, we can look at a population-time plot: on the X-axis we have time, and on the Y-axis we have the size of the risk set at a particular time point. Failure times associated to the event of interest can then be highlighted on the plot using red dots.

``` r
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

![](competingRisk_files/figure-markdown_github/poptime1-1.png)

We can right away draw a few conclusions from this plot: first of all, we get a sense of how quickly the size of the risk set changes over time. We also see that the incidence density is non-constant: most relapses occur before 15 months. Finally, we also see that the risk set keeps shrinking after the last event has occured; this could be due to either censoring or the competing event.

To get an idea of whether only relapse is responsible for the shrinking of the risk set in the first few months of follow-up, we can also keep track of how many events have occured at each time point:

``` r
# Plot at-risk and events
plot(0, type='n', xlim=c(0, max(ftime)), ylim=c(0, nobs), 
     xlab='Follow-up time', ylab='Population')
polygon(x = c(0,ftime[ord], max(ftime), 0), 
        y = c(0, yCoords[,2], 0, 0), 
        col = "firebrick3")
polygon(x = c(0, ftime[ord], ftime[rev(ord)], 0, 0),
        y = c(0, yCoords[,2], rev(yCoords[,2] + yCoords[,4]), nobs, 0), 
        col = "grey90")

# randomly move the cases vertically
moved_cases <- yCoords[cases[ord], 2] + yCoords[cases[ord], 4] * runif(sum(cases))
points((ftime[ord])[cases[ord]], moved_cases, pch=20, col="red", cex=1)
legend("topright", legend=c("Relapse", "At-risk"), 
       col=c("firebrick3", "grey90"),
       pch=15)
```

![](competingRisk_files/figure-markdown_github/poptime2-1.png)

Therefore, there is also censoring and loss due to competing events happening in the first few months. However, with this plot, we can't differentiate bwetween the two contributions. For this reason we can also keep track of the number of competing events at each time point:

``` r
plot(0, type='n', xlim=c(0, max(ftime)), ylim=c(0, nobs), 
     xlab='Follow-up time', ylab='Population')
polygon(x = c(0, max(ftime), max(ftime), 0),
        y = c(0, 0, nobs, nobs), col = "white")
# Event of interest
polygon(x = c(0,ftime[ord], max(ftime), 0), 
        y = c(0, yCoords[,2], 0, 0), 
        col = "firebrick3")
# Risk set
polygon(x = c(0, ftime[ord], ftime[rev(ord)], 0, 0),
        y = c(0, yCoords[,2], rev(yCoords[,2] + yCoords[,4]), nobs, 0), 
        col = "grey90")
# Competing event
polygon(x = c(0, ftime[ord], max(ftime), 0), 
        y = c(nobs, nobs - yCoords[,1], nobs, nobs), 
        col = "dodgerblue2")

# randomly move the cases vertically
moved_cases <- yCoords[cases[ord], 2] + yCoords[cases[ord], 4] * runif(sum(cases))
points((ftime[ord])[cases[ord]], moved_cases, pch=20, col="red", cex=1)
legend("topright", legend=c("Relapse", "Competing event", "At-risk"), 
       col=c("firebrick3", "dodgerblue2", "grey90"),
       pch=15)
```

![](competingRisk_files/figure-markdown_github/poptime3-1.png)

From this last plot, we can see that there is no censoring during the first 10 months. Moreover, we see that the last competing event occurs around 20 months. Putting all this information together, we have evidence of two types of patients: very sick patients who either relapse or have a competing event early on, and healthier patients who are eventually lost to follow-up.

Analysis
--------

We now turn to the analysis of this dataset. The population-time plots above give evidence of non-constant hazard; therefore, we will explicitely include time in the model. Note that we also include all other variables as possible confounders. First, we include time as a linear term:

``` r
library(casebase)
model1 <- fitSmoothHazard(Status ~ ftime + Sex + D + Phase + Source + Age, 
                          data = DT, 
                          ratio = 1000, 
                          type = "uniform", 
                          time = "ftime")
summary(model1)
```

    ## 
    ## Call:
    ## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
    ## 
    ## Pearson residuals:
    ##                         Min       1Q   Median        3Q    Max
    ## log(mu[,2]/mu[,1]) -0.07252 -0.02214 -0.01238 -0.004705 144.48
    ## log(mu[,3]/mu[,1]) -0.10874 -0.02844 -0.01188 -0.002703  63.81
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1  -3.471431   0.677926  -5.121 3.04e-07 ***
    ## (Intercept):2  -2.553408   0.456314  -5.596 2.20e-08 ***
    ## ftime:1        -0.069339   0.014637  -4.737 2.17e-06 ***
    ## ftime:2        -0.103005   0.018040  -5.710 1.13e-08 ***
    ## SexM:1         -0.283357   0.281482  -1.007 0.314097    
    ## SexM:2         -0.393539   0.234499  -1.678 0.093306 .  
    ## DAML:1         -0.615106   0.299529  -2.054 0.040016 *  
    ## DAML:2         -0.128529   0.275187  -0.467 0.640457    
    ## PhaseCR2:1      0.164757   0.465452   0.354 0.723360    
    ## PhaseCR2:2      0.284447   0.328909   0.865 0.387138    
    ## PhaseCR3:1      0.497277   0.689157   0.722 0.470557    
    ## PhaseCR3:2      0.247308   0.522180   0.474 0.635780    
    ## PhaseRelapse:1  1.425731   0.390068   3.655 0.000257 ***
    ## PhaseRelapse:2  0.760170   0.305914   2.485 0.012958 *  
    ## SourcePB:1      0.458579   0.569546   0.805 0.420725    
    ## SourcePB:2     -1.085761   0.353344  -3.073 0.002120 ** 
    ## Age:1          -0.006502   0.011923  -0.545 0.585557    
    ## Age:2           0.028097   0.009939   2.827 0.004699 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Dispersion Parameter for multinomial family:   1
    ## 
    ## Residual deviance: 2010.853 on 262244 degrees of freedom
    ## 
    ## Log-likelihood: -1005.426 on 262244 degrees of freedom
    ## 
    ## Number of iterations: 13 
    ## 
    ## Reference group is level  1  of the response

Because of the results in Turgeon *et al* \[-@turgeonCompRisk\], the standard errors we obtain from the multinomial logit fit are asymptotically correct, and therefore can be used to construct asymptotic confidence intervals.

From this summary, we see that time is indeed significant, as is Phase (only relapse vs. CR1). Interestingly, we see that the type of disease is only significant for the event of interest, whereas the type of transplant and the age of the patient are only significant for the competing event.

Next, we include the logarithm of time in the model (which leads to a Weibull hazard):

``` r
model2 <- fitSmoothHazard(Status ~ log(ftime) + Sex + D + Phase + Source + Age, 
                          data = DT, 
                          ratio = 1000, 
                          type = "uniform", 
                          time = "ftime")
summary(model2)
```

    ## 
    ## Call:
    ## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
    ## 
    ## Pearson residuals:
    ##                        Min       1Q   Median      3Q   Max
    ## log(mu[,2]/mu[,1]) -0.1586 -0.02172 -0.01488 -0.0112 91.50
    ## log(mu[,3]/mu[,1]) -0.2549 -0.02465 -0.01753 -0.0144 68.98
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1  -3.973373   0.697541  -5.696 1.22e-08 ***
    ## (Intercept):2  -3.072860   0.458942  -6.696 2.15e-11 ***
    ## log(ftime):1   -0.326885   0.068810  -4.751 2.03e-06 ***
    ## log(ftime):2   -0.409365   0.055783  -7.339 2.16e-13 ***
    ## SexM:1         -0.432316   0.291868  -1.481 0.138552    
    ## SexM:2         -0.508326   0.239301  -2.124 0.033653 *  
    ## DAML:1         -0.695149   0.301673  -2.304 0.021205 *  
    ## DAML:2         -0.159079   0.284358  -0.559 0.575868    
    ## PhaseCR2:1      0.256405   0.467125   0.549 0.583075    
    ## PhaseCR2:2      0.394798   0.329179   1.199 0.230396    
    ## PhaseCR3:1      0.437725   0.705766   0.620 0.535118    
    ## PhaseCR3:2      0.109218   0.530056   0.206 0.836752    
    ## PhaseRelapse:1  1.485346   0.392256   3.787 0.000153 ***
    ## PhaseRelapse:2  0.851896   0.307304   2.772 0.005569 ** 
    ## SourcePB:1      0.679091   0.594402   1.142 0.253256    
    ## SourcePB:2     -0.954314   0.364461  -2.618 0.008834 ** 
    ## Age:1          -0.003986   0.011672  -0.341 0.732740    
    ## Age:2           0.028748   0.009829   2.925 0.003445 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Dispersion Parameter for multinomial family:   1
    ## 
    ## Residual deviance: 2106.207 on 262244 degrees of freedom
    ## 
    ## Log-likelihood: -1053.103 on 262244 degrees of freedom
    ## 
    ## Number of iterations: 11 
    ## 
    ## Reference group is level  1  of the response

As we can see, the results are similar to the ones with a Gompertz hazard, although Sex is now significant for the competing event.

Finally, using splines, we can be quite flexible about the way the hazard depends on time:

``` r
model3 <- fitSmoothHazard(
    Status ~ splines::bs(ftime) + Sex + D + Phase + Source + Age, 
    data = DT, 
    ratio = 1000, 
    type = "uniform", 
    time = "ftime")
summary(model3)
```

    ## 
    ## Call:
    ## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
    ## 
    ## Pearson residuals:
    ##                         Min       1Q    Median         3Q   Max
    ## log(mu[,2]/mu[,1]) -0.06550 -0.02251 -0.012716 -2.422e-03 179.5
    ## log(mu[,3]/mu[,1]) -0.08964 -0.03003 -0.004575 -2.351e-06 116.1
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1          -3.691355   0.708817  -5.208 1.91e-07 ***
    ## (Intercept):2          -3.222215   0.510159  -6.316 2.68e-10 ***
    ## splines::bs(ftime)1:1  -0.065812   2.262489  -0.029 0.976794    
    ## splines::bs(ftime)1:2   7.050859   3.651841   1.931 0.053512 .  
    ## splines::bs(ftime)2:1 -16.116417   8.171818  -1.972 0.048587 *  
    ## splines::bs(ftime)2:2 -77.611389  25.533712  -3.040 0.002369 ** 
    ## splines::bs(ftime)3:1  -2.492098  10.137201  -0.246 0.805809    
    ## splines::bs(ftime)3:2  -2.618538  22.372792  -0.117 0.906827    
    ## SexM:1                 -0.312318   0.282712  -1.105 0.269280    
    ## SexM:2                 -0.432579   0.235026  -1.841 0.065686 .  
    ## DAML:1                 -0.637589   0.298442  -2.136 0.032648 *  
    ## DAML:2                 -0.150842   0.272504  -0.554 0.579894    
    ## PhaseCR2:1              0.177523   0.465556   0.381 0.702971    
    ## PhaseCR2:2              0.295998   0.330082   0.897 0.369857    
    ## PhaseCR3:1              0.535253   0.694359   0.771 0.440790    
    ## PhaseCR3:2              0.263867   0.525446   0.502 0.615543    
    ## PhaseRelapse:1          1.492639   0.395031   3.779 0.000158 ***
    ## PhaseRelapse:2          0.891996   0.310906   2.869 0.004117 ** 
    ## SourcePB:1              0.387227   0.572715   0.676 0.498961    
    ## SourcePB:2             -1.161639   0.356260  -3.261 0.001112 ** 
    ## Age:1                  -0.004674   0.011952  -0.391 0.695744    
    ## Age:2                   0.030146   0.010007   3.012 0.002591 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Dispersion Parameter for multinomial family:   1
    ## 
    ## Residual deviance: 1995.088 on 262240 degrees of freedom
    ## 
    ## Log-likelihood: -997.5441 on 262240 degrees of freedom
    ## 
    ## Number of iterations: 18 
    ## 
    ## Reference group is level  1  of the response

Again, we see that the results are quite similar for this third model.

### Absolute risk

We now look at the 2-year risk of relapse:

``` r
linearRisk <- absoluteRisk(object = model1, time = 24, newdata = DT[1:10,])
logRisk <- absoluteRisk(object = model2, time = 24, newdata = DT[1:10,])
splineRisk <- absoluteRisk(object = model3, time = 24, newdata = DT[1:10,])
```

``` r
plot(linearRisk[,1], logRisk[,1],
     xlab="Linear", ylab = "Log/Spline", pch=19,
     xlim=c(0,1), ylim=c(0,1), col='red')
points(linearRisk[,1], splineRisk[,1],
       col = 'blue', pch=19)
abline(a=0, b=1, lty=2, lwd=2)
legend("topleft", legend=c("Log", "Spline"),
       pch=19, col=c("red", "blue"))
```

![](competingRisk_files/figure-markdown_github/absRiskPlot-1.png)

As we can see, Model 1 and Model 2 give different absolute risk predictions, but the linear and the spline model actually give very similar results. We can also estimate the mean absolute risk for the entire dataset:

``` r
# The first column corresponds to the event of interest
mean(linearRisk[,1])
```

    ## [1] 0.204129

``` r
mean(logRisk[,1])
```

    ## [1] 0.1506065

``` r
mean(splineRisk[,1])
```

    ## [1] 0.1136662

Session information
-------------------

    ## R version 3.3.1 (2016-06-21)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.10
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] casebase_0.1.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.9      knitr_1.15.1     magrittr_1.5     splines_3.3.1   
    ##  [5] munsell_0.4.3    lattice_0.20-33  colorspace_1.3-1 stringr_1.2.0   
    ##  [9] plyr_1.8.4       tools_3.3.1      grid_3.3.1       data.table_1.9.6
    ## [13] gtable_0.2.0     htmltools_0.3.5  survival_2.39-5  yaml_2.1.14     
    ## [17] lazyeval_0.2.0   rprojroot_1.2    digest_0.6.12    assertthat_0.1  
    ## [21] tibble_1.2       Matrix_1.2-6     ggplot2_2.2.0    VGAM_1.0-2      
    ## [25] evaluate_0.10    rmarkdown_1.3    stringi_1.1.2    scales_0.4.1    
    ## [29] backports_1.0.5  stats4_3.3.1     chron_2.3-47

References
----------
