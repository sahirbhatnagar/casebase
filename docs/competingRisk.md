Data
----

We will use the same data that was used in Scrucca *et al* \[-@scrucca2010regression\]. The data is available on the main author's [website](http://www.stat.unipg.it/luca/R/).

``` r
set.seed(12345)
library(casebase)
data(bmtcrr)
head(bmtcrr)
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
nobs <- nrow(bmtcrr)
ftime <- bmtcrr$ftime
ord <- order(ftime, decreasing = FALSE)

# We split the person-moments in four categories:
# 1) at-risk
# 2) main event
# 3) competing event
# 4) censored
yCoords <- cbind(cumsum(bmtcrr[ord, "Status"] == 2), 
                 cumsum(bmtcrr[ord, "Status"] == 1),
                 cumsum(bmtcrr[ord, "Status"] == 0))
yCoords <- cbind(yCoords, nobs - rowSums(yCoords))

# Plot only at-risk
plot(0, type = 'n', xlim = c(0, max(ftime)), ylim = c(0, nobs), 
     xlab = 'Follow-up time', ylab = 'Population')
polygon(c(0, 0, ftime[ord], max(ftime), 0),
        c(0, nobs, yCoords[,4], 0, 0), col = "grey90")
cases <- bmtcrr[, "Status"] == 1

# randomly move the cases vertically
moved_cases <- yCoords[cases[ord], 4] * runif(sum(cases))
points((ftime[ord])[cases[ord]], moved_cases, pch = 20, 
       col = "red", cex = 1)
```

![](competingRisk_files/figure-markdown_github/poptime1-1.png)

We can right away draw a few conclusions from this plot: first of all, we get a sense of how quickly the size of the risk set changes over time. We also see that the incidence density is non-constant: most relapses occur before 15 months. Finally, we also see that the risk set keeps shrinking after the last event has occured; this could be due to either censoring or the competing event.

To get an idea of whether only relapse is responsible for the shrinking of the risk set in the first few months of follow-up, we can also keep track of how many events have occured at each time point:

``` r
# Plot at-risk and events
plot(0, type = 'n', xlim = c(0, max(ftime)), ylim = c(0, nobs), 
     xlab = 'Follow-up time', ylab = 'Population')
polygon(x = c(0,ftime[ord], max(ftime), 0), 
        y = c(0, yCoords[,2], 0, 0), 
        col = "firebrick3")
polygon(x = c(0, ftime[ord], ftime[rev(ord)], 0, 0),
        y = c(0, yCoords[,2], rev(yCoords[,2] + yCoords[,4]), nobs, 0), 
        col = "grey90")

# randomly move the cases vertically
moved_cases <- yCoords[cases[ord], 2] + yCoords[cases[ord], 4] * runif(sum(cases))
points((ftime[ord])[cases[ord]], moved_cases, pch = 20,
       col = "red", cex = 1)
legend("topright", legend = c("Relapse", "At-risk"), 
       col = c("firebrick3", "grey90"),
       pch = 15)
```

![](competingRisk_files/figure-markdown_github/poptime2-1.png)

Therefore, there is also censoring and loss due to competing events happening in the first few months. However, with this plot, we can't differentiate bwetween the two contributions. For this reason we can also keep track of the number of competing events at each time point:

``` r
plot(0, type = 'n', xlim = c(0, max(ftime)), ylim = c(0, nobs), 
     xlab = 'Follow-up time', ylab = 'Population')
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
points((ftime[ord])[cases[ord]], moved_cases, pch = 20,
       col = "red", cex = 1)
legend("topright", legend = c("Relapse", "Competing event", "At-risk"), 
       col = c("firebrick3", "dodgerblue2", "grey90"),
       pch = 15)
```

![](competingRisk_files/figure-markdown_github/poptime3-1.png)

From this last plot, we can see that there is no censoring during the first 10 months. Moreover, we see that the last competing event occurs around 20 months. Putting all this information together, we have evidence of two types of patients: very sick patients who either relapse or have a competing event early on, and healthier patients who are eventually lost to follow-up.

Analysis
--------

We now turn to the analysis of this dataset. The population-time plots above give evidence of non-constant hazard; therefore, we will explicitely include time in the model. Note that we also include all other variables as possible confounders. First, we include time as a linear term:

``` r
model1 <- fitSmoothHazard(Status ~ ftime + Sex + D + Phase + Source + Age, 
                          data = bmtcrr, 
                          ratio = 1000,
                          time = "ftime")
summary(model1)
```

    ## 
    ## Call:
    ## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
    ## 
    ## Pearson residuals:
    ##                         Min       1Q   Median        3Q    Max
    ## log(mu[,2]/mu[,1]) -0.07247 -0.02213 -0.01240 -0.004682 143.86
    ## log(mu[,3]/mu[,1]) -0.10769 -0.02841 -0.01182 -0.002664  63.76
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1  -3.493702   0.680148  -5.137 2.80e-07 ***
    ## (Intercept):2  -2.574417   0.459100  -5.608 2.05e-08 ***
    ## ftime:1        -0.069566   0.014665  -4.744 2.10e-06 ***
    ## ftime:2        -0.103680   0.018131  -5.718 1.08e-08 ***
    ## SexM:1         -0.279635   0.280845  -0.996 0.319401    
    ## SexM:2         -0.398319   0.234326  -1.700 0.089160 .  
    ## DAML:1         -0.612589   0.299105  -2.048 0.040553 *  
    ## DAML:2         -0.130847   0.273363  -0.479 0.632183    
    ## PhaseCR2:1      0.181035   0.465819   0.389 0.697544    
    ## PhaseCR2:2      0.301255   0.329598   0.914 0.360713    
    ## PhaseCR3:1      0.494073   0.689858   0.716 0.473871    
    ## PhaseCR3:2      0.224703   0.522988   0.430 0.667448    
    ## PhaseRelapse:1  1.434574   0.390435   3.674 0.000239 ***
    ## PhaseRelapse:2  0.768495   0.306073   2.511 0.012045 *  
    ## SourcePB:1      0.470195   0.567590   0.828 0.407441    
    ## SourcePB:2     -1.054388   0.350936  -3.005 0.002660 ** 
    ## Age:1          -0.006324   0.011931  -0.530 0.596064    
    ## Age:2           0.028052   0.009885   2.838 0.004542 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Dispersion Parameter for multinomial family:   1
    ## 
    ## Residual deviance: 2010.111 on 262244 degrees of freedom
    ## 
    ## Log-likelihood: -1005.056 on 262244 degrees of freedom
    ## 
    ## Number of iterations: 13 
    ## 
    ## Reference group is level  1  of the response

Because of the results in Turgeon *et al* \[-@turgeonCompRisk\], the standard errors we obtain from the multinomial logit fit are asymptotically correct, and therefore can be used to construct asymptotic confidence intervals.

From this summary, we see that time is indeed significant, as is Phase (only relapse vs. CR1). Interestingly, we see that the type of disease is only significant for the event of interest, whereas the type of transplant and the age of the patient are only significant for the competing event.

Next, we include the logarithm of time in the model (which leads to a Weibull hazard):

``` r
model2 <- fitSmoothHazard(Status ~ log(ftime) + Sex + D + Phase + Source + Age, 
                          data = bmtcrr, 
                          ratio = 1000, 
                          time = "ftime")
summary(model2)
```

    ## 
    ## Call:
    ## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
    ## 
    ## Pearson residuals:
    ##                        Min       1Q   Median       3Q   Max
    ## log(mu[,2]/mu[,1]) -0.1586 -0.02170 -0.01491 -0.01130 91.41
    ## log(mu[,3]/mu[,1]) -0.2765 -0.02455 -0.01753 -0.01442 68.65
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1  -3.976780   0.701316  -5.670 1.42e-08 ***
    ## (Intercept):2  -3.076452   0.460513  -6.680 2.38e-11 ***
    ## log(ftime):1   -0.331403   0.068786  -4.818 1.45e-06 ***
    ## log(ftime):2   -0.410863   0.055669  -7.380 1.58e-13 ***
    ## SexM:1         -0.406499   0.292048  -1.392 0.163956    
    ## SexM:2         -0.485123   0.239597  -2.025 0.042894 *  
    ## DAML:1         -0.674288   0.301859  -2.234 0.025497 *  
    ## DAML:2         -0.154872   0.284347  -0.545 0.585989    
    ## PhaseCR2:1      0.232940   0.467074   0.499 0.617976    
    ## PhaseCR2:2      0.357673   0.329746   1.085 0.278058    
    ## PhaseCR3:1      0.461212   0.711894   0.648 0.517072    
    ## PhaseCR3:2      0.118831   0.533390   0.223 0.823704    
    ## PhaseRelapse:1  1.462762   0.392443   3.727 0.000194 ***
    ## PhaseRelapse:2  0.850887   0.307694   2.765 0.005686 ** 
    ## SourcePB:1      0.671664   0.601670   1.116 0.264279    
    ## SourcePB:2     -0.968265   0.366793  -2.640 0.008295 ** 
    ## Age:1          -0.003776   0.011733  -0.322 0.747596    
    ## Age:2           0.029154   0.009940   2.933 0.003357 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Dispersion Parameter for multinomial family:   1
    ## 
    ## Residual deviance: 2106.329 on 262244 degrees of freedom
    ## 
    ## Log-likelihood: -1053.164 on 262244 degrees of freedom
    ## 
    ## Number of iterations: 11 
    ## 
    ## Reference group is level  1  of the response

As we can see, the results are similar to the ones with a Gompertz hazard, although Sex is now significant for the competing event.

Finally, using splines, we can be quite flexible about the way the hazard depends on time:

``` r
model3 <- fitSmoothHazard(
    Status ~ splines::bs(ftime) + Sex + D + Phase + Source + Age, 
    data = bmtcrr, 
    ratio = 1000, 
    time = "ftime")
summary(model3)
```

    ## 
    ## Call:
    ## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
    ## 
    ## Pearson residuals:
    ##                         Min       1Q    Median         3Q   Max
    ## log(mu[,2]/mu[,1]) -0.06497 -0.02241 -0.012772 -2.402e-03 178.2
    ## log(mu[,3]/mu[,1]) -0.08870 -0.03005 -0.004613 -2.360e-06 115.9
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1          -3.730220   0.707340  -5.274 1.34e-07 ***
    ## (Intercept):2          -3.260420   0.507168  -6.429 1.29e-10 ***
    ## splines::bs(ftime)1:1   0.049603   2.268962   0.022 0.982559    
    ## splines::bs(ftime)1:2   7.180721   3.646293   1.969 0.048916 *  
    ## splines::bs(ftime)2:1 -16.297428   8.159884  -1.997 0.045797 *  
    ## splines::bs(ftime)2:2 -77.667331  25.415907  -3.056 0.002244 ** 
    ## splines::bs(ftime)3:1  -2.511822  10.083023  -0.249 0.803273    
    ## splines::bs(ftime)3:2  -2.635626  22.561450  -0.117 0.907003    
    ## SexM:1                 -0.308727   0.282306  -1.094 0.274135    
    ## SexM:2                 -0.434290   0.234861  -1.849 0.064438 .  
    ## DAML:1                 -0.611650   0.299777  -2.040 0.041315 *  
    ## DAML:2                 -0.133471   0.273911  -0.487 0.626060    
    ## PhaseCR2:1              0.163280   0.465762   0.351 0.725914    
    ## PhaseCR2:2              0.295375   0.329884   0.895 0.370580    
    ## PhaseCR3:1              0.548603   0.694576   0.790 0.429622    
    ## PhaseCR3:2              0.292396   0.525374   0.557 0.577836    
    ## PhaseRelapse:1          1.490165   0.394595   3.776 0.000159 ***
    ## PhaseRelapse:2          0.897693   0.310735   2.889 0.003866 ** 
    ## SourcePB:1              0.416135   0.572748   0.727 0.467497    
    ## SourcePB:2             -1.152649   0.357278  -3.226 0.001254 ** 
    ## Age:1                  -0.005113   0.012054  -0.424 0.671442    
    ## Age:2                   0.030172   0.010050   3.002 0.002681 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Dispersion Parameter for multinomial family:   1
    ## 
    ## Residual deviance: 1995.49 on 262240 degrees of freedom
    ## 
    ## Log-likelihood: -997.745 on 262240 degrees of freedom
    ## 
    ## Number of iterations: 18 
    ## 
    ## Reference group is level  1  of the response

Again, we see that the results are quite similar for this third model.

### Absolute risk

We now look at the 2-year risk of relapse:

``` r
linearRisk <- absoluteRisk(object = model1, time = 24, newdata = bmtcrr[1:10,])
logRisk <- absoluteRisk(object = model2, time = 24, newdata = bmtcrr[1:10,])
splineRisk <- absoluteRisk(object = model3, time = 24, newdata = bmtcrr[1:10,])
```

``` r
plot(linearRisk[,1], logRisk[,1],
     xlab = "Linear", ylab = "Log/Spline", pch = 19,
     xlim = c(0,1), ylim = c(0,1), col = 'red')
points(linearRisk[,1], splineRisk[,1],
       col = 'blue', pch = 19)
abline(a = 0, b = 1, lty = 2, lwd = 2)
legend("topleft", legend = c("Log", "Spline"),
       pch = 19, col = c("red", "blue"))
```

![](competingRisk_files/figure-markdown_github/absRiskPlot-1.png)

As we can see, Model 1 and Model 2 give different absolute risk predictions, but the linear and the spline model actually give very similar results. We can also estimate the mean absolute risk for the entire dataset:

``` r
# The first column corresponds to the event of interest
mean(linearRisk[,1])
```

    ## [1] 0.1422626

``` r
mean(logRisk[,1])
```

    ## [1] 0.1816989

``` r
mean(splineRisk[,1])
```

    ## [1] 0.1393753

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
