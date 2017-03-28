Data
----

We will use the same data that was used in Scrucca *et al* \[-@scrucca2010regression\]. The data is available on the main author's [website](http://www.stat.unipg.it/luca/R/).

``` r
set.seed(12345)

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
    ##                         Min       1Q   Median       3Q   Max
    ## log(mu[,2]/mu[,1]) -0.03511 -0.02312 -0.01946 -0.01588 84.80
    ## log(mu[,3]/mu[,1]) -0.03693 -0.02680 -0.02337 -0.02017 57.29
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1  -4.7086629  0.6860839  -6.863 6.74e-12 ***
    ## (Intercept):2  -3.7753463  0.4712041  -8.012 1.13e-15 ***
    ## ftime:1        -0.0074411  0.0095213  -0.782   0.4345    
    ## ftime:2        -0.0234739  0.0118882  -1.975   0.0483 *  
    ## SexM:1          0.1173418  0.2766552   0.424   0.6715    
    ## SexM:2         -0.1317348  0.2355512  -0.559   0.5760    
    ## DAML:1         -0.2615793  0.3039472  -0.861   0.3895    
    ## DAML:2          0.1176449  0.2671575   0.440   0.6597    
    ## PhaseCR2:1      0.1072881  0.4623358   0.232   0.8165    
    ## PhaseCR2:2      0.1036405  0.3310731   0.313   0.7542    
    ## PhaseCR3:1      0.2796638  0.6741410   0.415   0.6783    
    ## PhaseCR3:2      0.1126718  0.5154568   0.219   0.8270    
    ## PhaseRelapse:1  0.8054620  0.3804372   2.117   0.0342 *  
    ## PhaseRelapse:2 -0.0007109  0.3000402  -0.002   0.9981    
    ## SourcePB:1      0.7456485  0.5311726   1.404   0.1604    
    ## SourcePB:2     -0.6471805  0.3288367  -1.968   0.0491 *  
    ## Age:1          -0.0151611  0.0118020  -1.285   0.1989    
    ## Age:2           0.0186593  0.0098371   1.897   0.0579 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Dispersion Parameter for multinomial family:   1
    ## 
    ## Residual deviance: 2221.211 on 262244 degrees of freedom
    ## 
    ## Log-likelihood: -1110.605 on 262244 degrees of freedom
    ## 
    ## Number of iterations: 11 
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
    ##                         Min       1Q   Median       3Q   Max
    ## log(mu[,2]/mu[,1]) -0.05154 -0.02288 -0.01882 -0.01490 79.90
    ## log(mu[,3]/mu[,1]) -0.04174 -0.02653 -0.02323 -0.01983 64.38
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1  -5.169337   0.689581  -7.496 6.56e-14 ***
    ## (Intercept):2  -4.150096   0.469486  -8.840  < 2e-16 ***
    ## log(ftime):1    0.275468   0.087632   3.143  0.00167 ** 
    ## log(ftime):2    0.136730   0.074470   1.836  0.06635 .  
    ## SexM:1         -0.101932   0.282880  -0.360  0.71860    
    ## SexM:2         -0.277052   0.235713  -1.175  0.23984    
    ## DAML:1         -0.462982   0.304138  -1.522  0.12794    
    ## DAML:2          0.022211   0.273592   0.081  0.93530    
    ## PhaseCR2:1      0.135328   0.463264   0.292  0.77020    
    ## PhaseCR2:2      0.183499   0.329272   0.557  0.57733    
    ## PhaseCR3:1      0.350681   0.680687   0.515  0.60642    
    ## PhaseCR3:2      0.137675   0.516486   0.267  0.78981    
    ## PhaseRelapse:1  1.106652   0.389873   2.838  0.00453 ** 
    ## PhaseRelapse:2  0.304955   0.309689   0.985  0.32477    
    ## SourcePB:1      0.705531   0.549936   1.283  0.19952    
    ## SourcePB:2     -0.737839   0.339630  -2.172  0.02982 *  
    ## Age:1          -0.010259   0.011751  -0.873  0.38264    
    ## Age:2           0.021852   0.009808   2.228  0.02587 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Dispersion Parameter for multinomial family:   1
    ## 
    ## Residual deviance: 2214.089 on 262244 degrees of freedom
    ## 
    ## Log-likelihood: -1107.044 on 262244 degrees of freedom
    ## 
    ## Number of iterations: 10 
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
    ##                         Min       1Q   Median       3Q   Max
    ## log(mu[,2]/mu[,1]) -0.05787 -0.02362 -0.01828 -0.01349 103.9
    ## log(mu[,3]/mu[,1]) -0.05795 -0.02821 -0.02077 -0.01506  93.3
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1          -5.322353   0.701356  -7.589 3.23e-14 ***
    ## (Intercept):2          -4.806549   0.505056  -9.517  < 2e-16 ***
    ## splines::bs(ftime)1:1  10.051640   2.514058   3.998 6.38e-05 ***
    ## splines::bs(ftime)1:2  18.542863   3.706641   5.003 5.66e-07 ***
    ## splines::bs(ftime)2:1 -27.155285   9.621508  -2.822 0.004767 ** 
    ## splines::bs(ftime)2:2 -99.768340  26.045064  -3.831 0.000128 ***
    ## splines::bs(ftime)3:1   2.032730   7.508245   0.271 0.786596    
    ## splines::bs(ftime)3:2   0.716355  21.435757   0.033 0.973341    
    ## SexM:1                  0.010810   0.277099   0.039 0.968881    
    ## SexM:2                 -0.224871   0.236840  -0.949 0.342385    
    ## DAML:1                 -0.414197   0.302615  -1.369 0.171085    
    ## DAML:2                 -0.005705   0.265368  -0.021 0.982847    
    ## PhaseCR2:1              0.040206   0.463103   0.087 0.930815    
    ## PhaseCR2:2              0.121606   0.330203   0.368 0.712666    
    ## PhaseCR3:1              0.535662   0.679935   0.788 0.430805    
    ## PhaseCR3:2              0.274984   0.518349   0.531 0.595765    
    ## PhaseRelapse:1          1.086171   0.391777   2.772 0.005564 ** 
    ## PhaseRelapse:2          0.328237   0.310983   1.055 0.291205    
    ## SourcePB:1              0.507192   0.537842   0.943 0.345674    
    ## SourcePB:2             -0.892167   0.336761  -2.649 0.008067 ** 
    ## Age:1                  -0.012264   0.012064  -1.017 0.309381    
    ## Age:2                   0.022048   0.009932   2.220 0.026427 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Dispersion Parameter for multinomial family:   1
    ## 
    ## Residual deviance: 2157.524 on 262240 degrees of freedom
    ## 
    ## Log-likelihood: -1078.762 on 262240 degrees of freedom
    ## 
    ## Number of iterations: 16 
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
