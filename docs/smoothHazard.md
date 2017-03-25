Methodological details
----------------------

Case-base sampling was proposed by Hanley & Miettinen \[-@hanley2009fitting\] as a way to fit smooth-in-time parametric hazard functions via logistic regression. The main idea, which was first proposed by Mantel \[-@mantel1973synthetic\] and then later developped by Efron \[-@efron1977efficiency\], is to sample person-moments, i.e. discrete time points along an subject's follow-up time, in order to construct a base series against which the case series can be compared.

This approach allows the explicit inclusion of the time variable into the model, which enables the user to fit a wide class of parametric hazard functions. For example, including time linearly recovers the Gompertz hazard, whereas including time *logarithmically* recovers the Weibull hazard; not including time at all corresponds to the exponential hazard.

The theoretical properties of this approach have been studied in Saarela & Arjas \[-@saarela2015non\] and Saarela \[-@saarela2015case\].

First example
-------------

The first example we discuss uses the well-known `veteran` dataset, which is part of the `survival` package. As we can see below, there is almost no censoring, and therefore we can get a good visual representation of the survival function:

``` r
library(survival)
data(veteran)
table(veteran$status)
```

    ## 
    ##   0   1 
    ##   9 128

``` r
evtimes <- veteran$time[veteran$status == 1]
hist(evtimes, nclass = 30, main = '', xlab = 'Survival time (days)', 
     col = 'gray90', probability = TRUE)
tgrid <- seq(0, 1000, by = 10)
lines(tgrid, dexp(tgrid, rate = 1.0/mean(evtimes)), 
      lwd = 2, lty = 2, col = 'red')
```

![](smoothHazard_files/figure-markdown_github/unnamed-chunk-1-1.png)

As we can see, the empirical survival function ressembles an exponential distribution.

We will first try to estimate the hazard function parametrically using some well-known regression routines. But first, we will reformat the data slightly.

``` r
veteran$prior <- factor(veteran$prior, levels = c(0, 10))
veteran$celltype <- factor(veteran$celltype, 
                           levels = c('large', 'squamous', 'smallcell', 'adeno'))
veteran$trt <- factor(veteran$trt, levels = c(1, 2))
```

Using the `eha` package, we can fit a Weibull form, with different values of the shape parameter. For `shape = 1`, we get an exponential distribution:

``` r
pacman::p_load(eha)
y <- with(veteran, Surv(time, status))

model1 <- weibreg(y ~ karno + diagtime + age + prior + celltype + trt, 
                  data = veteran, shape = 1)
summary(model1)
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

If we take `shape = 0`, the shape parameter is estimated along with the regression coefficients:

``` r
model2 <- weibreg(y ~ karno + diagtime + age + prior + celltype + trt, 
                  data = veteran, shape = 0)
summary(model2)
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

Finally, we can also fit a Cox proportional hazard:

``` r
model3 <- coxph(y ~ karno + diagtime + age + prior + celltype + trt, 
                data = veteran)
summary(model3)
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

As we can see, all three models are significant, and they give similar information: `karno` and `celltype` are significant predictors, both treatment is not.

The method available in this package makes use of *case-base sampling*. That is, person-moments are randomly sampled across the entire follow-up time, with some moments corresponding to cases and others to controls. By sampling person-moments instead of individuals, we can then use logistic regression to fit smooth-in-time parametric hazard functions. See the previous section for more details.

First, we will look at the follow-up time by using population-time plots:

``` r
nobs <- nrow(y)
ftime <- veteran$time
ord <- order(ftime, decreasing = TRUE)
plot(0, type = 'n', xlim = c(0, max(ftime)), ylim = c(0, nobs), 
     xlab = 'Follow-up time', ylab = 'Population')
segments(rep(0.0, nobs), 1:nobs, ftime[ord], 1:nobs, col = 'gray25')
cases <- veteran$status == 1
points((ftime[ord])[cases[ord]], (1:nobs)[cases[ord]], pch = 20, col = 'red', cex = 0.5)
```

![](smoothHazard_files/figure-markdown_github/unnamed-chunk-6-1.png)

Population-time plots are a useful way of visualizing the total follow-up experience, where individuals appear on the y-axis, and follow-up time on the x-axis; each individual's follow-up time is represented by a gray line segment. For convenience, we have ordered the patients according to their time-to-event, and each event is represented by a red dot. The censored observations (of which there is only a few) correspond to the grey lines which do not end with a red dot.

Next, we use case-base sampling to fit a parametric hazard function via logistic regression. First, we will include time as a linear term; as noted above, this corresponds to an Gompertz hazard.

``` r
library(casebase)
model4 <- fitSmoothHazard(status ~ time + karno + diagtime + age + prior +
             celltype + trt, data = veteran, ratio = 100, type = "uniform")
```

    ## 'time' will be used as the time variable

``` r
summary(model4)
```

    ## 
    ## Call:
    ## glm(formula = formula, family = binomial, data = sampleData)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -0.4076  -0.1506  -0.1205  -0.1012   3.3914  
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -2.8090170  0.7203905  -3.899 9.65e-05 ***
    ## time               0.0003933  0.0006454   0.609   0.5422    
    ## karno             -0.0302912  0.0052240  -5.798 6.69e-09 ***
    ## diagtime           0.0023264  0.0091108   0.255   0.7985    
    ## age               -0.0071628  0.0092697  -0.773   0.4397    
    ## prior10           -0.0393741  0.2286357  -0.172   0.8633    
    ## celltypesquamous  -0.4025799  0.2843369  -1.416   0.1568    
    ## celltypesmallcell  0.4656078  0.2626241   1.773   0.0762 .  
    ## celltypeadeno      0.7408300  0.3025105   2.449   0.0143 *  
    ## trt2               0.2133287  0.2018410   1.057   0.2906    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1436.2  on 12927  degrees of freedom
    ## Residual deviance: 1367.1  on 12918  degrees of freedom
    ## AIC: 1387.1
    ## 
    ## Number of Fisher Scoring iterations: 8

Since the output object from `fitSmoothHazard` inherits from the `glm` class, we see a familiar result when using the function `summary`.

The main purpose of fitting smooth hazard functions is that it is then relatively easy to compute absolute risks. For example, we can use the function `absoluteRisk` to compute the mean absolute risk at 90 days, which can then be compared to the empirical measure.

``` r
absoluteRisk(object = model4, time = 90)
```

    ## [1] 0.5758602

``` r
mean(ftime <= 90)
```

    ## [1] 0.5547445

We can also fit a Weibull hazard by using a logarithmic term for time:

``` r
model5 <- fitSmoothHazard(status ~ log(time) + karno + diagtime + age + prior +
             celltype + trt, data = veteran, ratio = 100, type = "uniform")
```

    ## 'time' will be used as the time variable

``` r
summary(model5)
```

    ## 
    ## Call:
    ## glm(formula = formula, family = binomial, data = sampleData)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -0.4374  -0.1525  -0.1211  -0.0974   3.4403  
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -2.869604   0.750754  -3.822 0.000132 ***
    ## log(time)          0.062363   0.072533   0.860 0.389908    
    ## karno             -0.032509   0.005475  -5.938 2.89e-09 ***
    ## diagtime          -0.001032   0.009346  -0.110 0.912091    
    ## age               -0.007281   0.009258  -0.786 0.431599    
    ## prior10            0.089708   0.229979   0.390 0.696483    
    ## celltypesquamous  -0.476523   0.279066  -1.708 0.087718 .  
    ## celltypesmallcell  0.432164   0.264012   1.637 0.101650    
    ## celltypeadeno      0.736092   0.302782   2.431 0.015053 *  
    ## trt2               0.250236   0.203763   1.228 0.219420    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1436.2  on 12927  degrees of freedom
    ## Residual deviance: 1365.9  on 12918  degrees of freedom
    ## AIC: 1385.9
    ## 
    ## Number of Fisher Scoring iterations: 8

With case-base sampling, it is straightforward to fit a semi-parametric hazard function using splines, which can then be used to estimate the mean absolute risk.

``` r
# Fit a spline for time
library(splines)
model6 <- fitSmoothHazard(status ~ bs(time) + karno + diagtime + age + prior +
             celltype + trt, data = veteran, ratio = 100, type = "uniform")
```

    ## 'time' will be used as the time variable

``` r
summary(model6)
```

    ## 
    ## Call:
    ## glm(formula = formula, family = binomial, data = sampleData)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -0.4636  -0.1516  -0.1180  -0.0944   3.5296  
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -2.843120   0.725929  -3.917 8.98e-05 ***
    ## bs(time)1          1.617602   1.021959   1.583  0.11346    
    ## bs(time)2         -2.553143   1.736336  -1.470  0.14145    
    ## bs(time)3          1.806882   1.016930   1.777  0.07560 .  
    ## karno             -0.033919   0.005392  -6.290 3.17e-10 ***
    ## diagtime           0.004330   0.009963   0.435  0.66381    
    ## age               -0.007121   0.009311  -0.765  0.44441    
    ## prior10           -0.018463   0.234631  -0.079  0.93728    
    ## celltypesquamous  -0.380942   0.283311  -1.345  0.17875    
    ## celltypesmallcell  0.484411   0.264278   1.833  0.06681 .  
    ## celltypeadeno      0.820671   0.302580   2.712  0.00668 ** 
    ## trt2               0.296401   0.206425   1.436  0.15104    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1436.2  on 12927  degrees of freedom
    ## Residual deviance: 1359.5  on 12916  degrees of freedom
    ## AIC: 1383.5
    ## 
    ## Number of Fisher Scoring iterations: 8

``` r
absoluteRisk(object = model6, time = 90)
```

    ## [1] 0.5755579

As we can see from the summary, there is little evidence that splines actually improve the fit. Moreover, we can see that estimated individual absolute risks are essentially the same when using either a linear term or splines:

``` r
linearRisk <- absoluteRisk(object = model4, time = 90, newdata = veteran)
splineRisk <- absoluteRisk(object = model6, time = 90, newdata = veteran)

plot(linearRisk, splineRisk,
     xlab = "Linear", ylab = "Splines", pch = 19)
abline(a = 0, b = 1, lty = 2, lwd = 2, col = 'red')
```

![](smoothHazard_files/figure-markdown_github/unnamed-chunk-11-1.png)

These last three models give similar information as the first three, i.e. the main predictors for the hazard are `karno` and `celltype`, with treatment being non-significant. Moreover, by explicitely including the time variable in the formula, we see that it is not significant; this is evidence that the true hazard is exponential.

Finally, we can look at the estimates of the coefficients for the Cox model, as well as the last three models (CB stands for "case-base"):

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Cox model
</th>
<th style="text-align:right;">
CB linear
</th>
<th style="text-align:right;">
CB log-linear
</th>
<th style="text-align:right;">
CB splines
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
karno
</td>
<td style="text-align:right;">
-0.0328
</td>
<td style="text-align:right;">
-0.0303
</td>
<td style="text-align:right;">
-0.0325
</td>
<td style="text-align:right;">
-0.0339
</td>
</tr>
<tr>
<td style="text-align:left;">
diagtime
</td>
<td style="text-align:right;">
0.0001
</td>
<td style="text-align:right;">
0.0023
</td>
<td style="text-align:right;">
-0.0010
</td>
<td style="text-align:right;">
0.0043
</td>
</tr>
<tr>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
-0.0087
</td>
<td style="text-align:right;">
-0.0072
</td>
<td style="text-align:right;">
-0.0073
</td>
<td style="text-align:right;">
-0.0071
</td>
</tr>
<tr>
<td style="text-align:left;">
prior10
</td>
<td style="text-align:right;">
0.0716
</td>
<td style="text-align:right;">
-0.0394
</td>
<td style="text-align:right;">
0.0897
</td>
<td style="text-align:right;">
-0.0185
</td>
</tr>
<tr>
<td style="text-align:left;">
celltypesquamous
</td>
<td style="text-align:right;">
-0.4013
</td>
<td style="text-align:right;">
-0.4026
</td>
<td style="text-align:right;">
-0.4765
</td>
<td style="text-align:right;">
-0.3809
</td>
</tr>
<tr>
<td style="text-align:left;">
celltypesmallcell
</td>
<td style="text-align:right;">
0.4603
</td>
<td style="text-align:right;">
0.4656
</td>
<td style="text-align:right;">
0.4322
</td>
<td style="text-align:right;">
0.4844
</td>
</tr>
<tr>
<td style="text-align:left;">
celltypeadeno
</td>
<td style="text-align:right;">
0.7948
</td>
<td style="text-align:right;">
0.7408
</td>
<td style="text-align:right;">
0.7361
</td>
<td style="text-align:right;">
0.8207
</td>
</tr>
<tr>
<td style="text-align:left;">
trt2
</td>
<td style="text-align:right;">
0.2946
</td>
<td style="text-align:right;">
0.2133
</td>
<td style="text-align:right;">
0.2502
</td>
<td style="text-align:right;">
0.2964
</td>
</tr>
</tbody>
</table>
Session information
-------------------

    ## R version 3.3.1 (2016-06-21)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.10
    ## 
    ## attached base packages:
    ## [1] splines   stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] casebase_0.1.0  eha_2.4-4       survival_2.39-5
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.9      knitr_1.15.1     magrittr_1.5     munsell_0.4.3   
    ##  [5] colorspace_1.3-1 lattice_0.20-33  highr_0.6        plyr_1.8.4      
    ##  [9] stringr_1.2.0    tools_3.3.1      grid_3.3.1       data.table_1.9.6
    ## [13] gtable_0.2.0     pacman_0.4.1     htmltools_0.3.5  assertthat_0.1  
    ## [17] lazyeval_0.2.0   yaml_2.1.14      rprojroot_1.2    digest_0.6.12   
    ## [21] tibble_1.2       Matrix_1.2-6     ggplot2_2.2.0    VGAM_1.0-2      
    ## [25] evaluate_0.10    rmarkdown_1.3    stringi_1.1.2    scales_0.4.1    
    ## [29] backports_1.0.5  stats4_3.3.1     chron_2.3-47

References
----------
