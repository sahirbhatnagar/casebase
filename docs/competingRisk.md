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
                          ratio = 100,
                          time = "ftime")
summary(model1)
```

    ## 
    ## Call:
    ## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
    ## 
    ## 
    ## Pearson residuals:
    ##                        Min       1Q   Median        3Q   Max
    ## log(mu[,2]/mu[,1]) -0.2240 -0.06967 -0.03872 -0.014616 46.24
    ## log(mu[,3]/mu[,1]) -0.3356 -0.08920 -0.03671 -0.008251 20.23
    ## 
    ## Coefficients: 
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1  -3.527146   0.685168  -5.148 2.63e-07 ***
    ## (Intercept):2  -2.648451   0.463012  -5.720 1.06e-08 ***
    ## ftime:1        -0.070927   0.014929  -4.751 2.02e-06 ***
    ## ftime:2        -0.105177   0.018349  -5.732 9.93e-09 ***
    ## SexM:1         -0.289067   0.283217  -1.021 0.307418    
    ## SexM:2         -0.382981   0.236935  -1.616 0.106008    
    ## DAML:1         -0.575749   0.299617  -1.922 0.054654 .  
    ## DAML:2         -0.100149   0.274099  -0.365 0.714833    
    ## PhaseCR2:1      0.186766   0.467042   0.400 0.689237    
    ## PhaseCR2:2      0.286425   0.332270   0.862 0.388673    
    ## PhaseCR3:1      0.586630   0.696521   0.842 0.399660    
    ## PhaseCR3:2      0.310781   0.530986   0.585 0.558353    
    ## PhaseRelapse:1  1.448907   0.391878   3.697 0.000218 ***
    ## PhaseRelapse:2  0.792938   0.307933   2.575 0.010023 *  
    ## SourcePB:1      0.456442   0.571108   0.799 0.424162    
    ## SourcePB:2     -1.013983   0.355666  -2.851 0.004359 ** 
    ## Age:1          -0.005242   0.011917  -0.440 0.660007    
    ## Age:2           0.028597   0.009929   2.880 0.003976 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Residual deviance: 1409.076 on 26444 degrees of freedom
    ## 
    ## Log-likelihood: -704.5378 on 26444 degrees of freedom
    ## 
    ## Number of iterations: 10 
    ## 
    ## Reference group is level  1  of the response

Because of the results in Turgeon *et al* \[-@turgeonCompRisk\], the standard errors we obtain from the multinomial logit fit are asymptotically correct, and therefore can be used to construct asymptotic confidence intervals.

From this summary, we see that time is indeed significant, as is Phase (only relapse vs. CR1). Interestingly, we see that the type of disease is only significant for the event of interest, whereas the type of transplant and the age of the patient are only significant for the competing event.

Next, we include the logarithm of time in the model (which leads to a Weibull hazard):

``` r
model2 <- fitSmoothHazard(Status ~ log(ftime) + Sex + D + Phase + Source + Age, 
                          data = bmtcrr, 
                          ratio = 100, 
                          time = "ftime")
summary(model2)
```

    ## 
    ## Call:
    ## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
    ## 
    ## 
    ## Pearson residuals:
    ##                        Min       1Q   Median       3Q   Max
    ## log(mu[,2]/mu[,1]) -0.3728 -0.06889 -0.04686 -0.03489 29.01
    ## log(mu[,3]/mu[,1]) -0.7720 -0.07784 -0.05517 -0.04490 22.03
    ## 
    ## Coefficients: 
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1  -3.947939   0.704805  -5.601 2.13e-08 ***
    ## (Intercept):2  -3.024184   0.464824  -6.506 7.71e-11 ***
    ## log(ftime):1   -0.341809   0.070171  -4.871 1.11e-06 ***
    ## log(ftime):2   -0.426477   0.057517  -7.415 1.22e-13 ***
    ## SexM:1         -0.443805   0.292859  -1.515 0.129666    
    ## SexM:2         -0.519251   0.241471  -2.150 0.031526 *  
    ## DAML:1         -0.740607   0.301246  -2.458 0.013953 *  
    ## DAML:2         -0.141243   0.287950  -0.491 0.623771    
    ## PhaseCR2:1      0.305084   0.469113   0.650 0.515471    
    ## PhaseCR2:2      0.443158   0.331173   1.338 0.180849    
    ## PhaseCR3:1      0.467723   0.707953   0.661 0.508824    
    ## PhaseCR3:2      0.117413   0.532943   0.220 0.825629    
    ## PhaseRelapse:1  1.472912   0.394988   3.729 0.000192 ***
    ## PhaseRelapse:2  0.834838   0.310771   2.686 0.007224 ** 
    ## SourcePB:1      0.696141   0.598801   1.163 0.245009    
    ## SourcePB:2     -0.971909   0.374933  -2.592 0.009536 ** 
    ## Age:1          -0.003405   0.011665  -0.292 0.770361    
    ## Age:2           0.028570   0.009802   2.915 0.003561 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Residual deviance: 1503.138 on 26444 degrees of freedom
    ## 
    ## Log-likelihood: -751.5688 on 26444 degrees of freedom
    ## 
    ## Number of iterations: 8 
    ## 
    ## Reference group is level  1  of the response

As we can see, the results are similar to the ones with a Gompertz hazard, although Sex is now significant for the competing event.

Finally, using splines, we can be quite flexible about the way the hazard depends on time:

``` r
model3 <- fitSmoothHazard(
    Status ~ splines::bs(ftime) + Sex + D + Phase + Source + Age, 
    data = bmtcrr, 
    ratio = 100, 
    time = "ftime")
summary(model3)
```

    ## 
    ## Call:
    ## vglm(formula = formula, family = multinomial(refLevel = 1), data = sampleData)
    ## 
    ## 
    ## Pearson residuals:
    ##                        Min       1Q   Median         3Q   Max
    ## log(mu[,2]/mu[,1]) -0.2112 -0.07072 -0.03924 -7.966e-03 55.79
    ## log(mu[,3]/mu[,1]) -0.2810 -0.09476 -0.01408 -9.698e-06 34.23
    ## 
    ## Coefficients: 
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept):1          -3.756937   0.707736  -5.308 1.11e-07 ***
    ## (Intercept):2          -3.230327   0.512122  -6.308 2.83e-10 ***
    ## splines::bs(ftime)1:1  -0.157959   2.277313  -0.069 0.944701    
    ## splines::bs(ftime)1:2   6.639515   3.695359   1.797 0.072381 .  
    ## splines::bs(ftime)2:1 -15.907763   8.133390  -1.956 0.050482 .  
    ## splines::bs(ftime)2:2 -75.091363  25.403026  -2.956 0.003117 ** 
    ## splines::bs(ftime)3:1  -2.415094  10.618108  -0.227 0.820073    
    ## splines::bs(ftime)3:2  -2.157054  22.830019  -0.094 0.924725    
    ## SexM:1                 -0.327386   0.283998  -1.153 0.249001    
    ## SexM:2                 -0.483204   0.237367  -2.036 0.041782 *  
    ## DAML:1                 -0.611398   0.301775  -2.026 0.042764 *  
    ## DAML:2                 -0.112716   0.274870  -0.410 0.681753    
    ## PhaseCR2:1              0.215153   0.468799   0.459 0.646273    
    ## PhaseCR2:2              0.418330   0.333261   1.255 0.209383    
    ## PhaseCR3:1              0.481796   0.694753   0.693 0.488010    
    ## PhaseCR3:2              0.214904   0.528079   0.407 0.684042    
    ## PhaseRelapse:1          1.519795   0.395642   3.841 0.000122 ***
    ## PhaseRelapse:2          0.941056   0.313721   3.000 0.002703 ** 
    ## SourcePB:1              0.474662   0.573470   0.828 0.407839    
    ## SourcePB:2             -1.163727   0.359718  -3.235 0.001216 ** 
    ## Age:1                  -0.005299   0.012045  -0.440 0.660024    
    ## Age:2                   0.030127   0.009992   3.015 0.002569 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Number of linear predictors:  2 
    ## 
    ## Names of linear predictors: log(mu[,2]/mu[,1]), log(mu[,3]/mu[,1])
    ## 
    ## Residual deviance: 1392.237 on 26440 degrees of freedom
    ## 
    ## Log-likelihood: -696.1185 on 26440 degrees of freedom
    ## 
    ## Number of iterations: 16 
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
plot(linearRisk[,,1], logRisk[,,1],
     xlab = "Linear", ylab = "Log/Spline", pch = 19,
     xlim = c(0,1), ylim = c(0,1), col = 'red')
points(linearRisk[,,1], splineRisk[,,1],
       col = 'blue', pch = 19)
abline(a = 0, b = 1, lty = 2, lwd = 2)
legend("topleft", legend = c("Log", "Spline"),
       pch = 19, col = c("red", "blue"))
```

![](competingRisk_files/figure-markdown_github/absRiskPlot-1.png)

We can also estimate the mean absolute risk for the entire dataset:

``` r
mean(linearRisk[,,1])
```

    ## [1] 0.2012265

``` r
mean(logRisk[,,1])
```

    ## [1] 0.1453269

``` r
mean(splineRisk[,,1])
```

    ## [1] 0.105314

Session information
-------------------

    ## R version 3.3.3 (2017-03-06)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.2 LTS
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] casebase_0.1.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.10       knitr_1.15.1       magrittr_1.5      
    ##  [4] splines_3.3.3      munsell_0.4.3      lattice_0.20-35   
    ##  [7] colorspace_1.3-2   stringr_1.2.0      plyr_1.8.4        
    ## [10] tools_3.3.3        grid_3.3.3         data.table_1.10.4 
    ## [13] gtable_0.2.0       htmltools_0.3.6    survival_2.40-1   
    ## [16] yaml_2.1.14        lazyeval_0.2.0     rprojroot_1.2     
    ## [19] digest_0.6.12      tibble_1.3.0       Matrix_1.2-8      
    ## [22] ggplot2_2.2.1      VGAM_1.0-3         evaluate_0.10     
    ## [25] rmarkdown_1.3.9003 stringi_1.1.5      scales_0.4.1      
    ## [28] backports_1.0.5    stats4_3.3.3
