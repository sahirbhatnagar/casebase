---
author: Maxime
date: 'January 22, 2017'
output:
  md_document:
    variant: markdown
title: data
---

Methodological details
----------------------

Case-base sampling was proposed by [Hanley and Miettinen,
2009](https://github.com/sahirbhatnagar/casebase/blob/master/references/Hanley_Miettinen-2009-Inter_J_of_Biostats.pdf)
as a way to fit smooth-in-time parametric hazard functions via logistic
regression. The main idea, which was first proposed by Mantel, 1973 and
then later developped by Efron, 1977, is to sample person-moments, i.e.
discrete time points along an subject's follow-up time, in order to
construct a base series against which the case series can be compared.

This approach allows the explicit inclusion of the time variable into
the model, which enables the user to fit a wide class of parametric
hazard functions. For example, including time linearly recovers the
Gompertz hazard, whereas including time *logarithmically* recovers the
Weibull hazard; not including time at all corresponds to the exponential
hazard.

The theoretical properties of this approach have been studied in
[Saarela and Arjas,
2015](https://github.com/sahirbhatnagar/casebase/blob/master/references/Saarela_et_al-2015-Scandinavian_Journal_of_Statistics.pdf)
and [Saarela,
2015](https://github.com/sahirbhatnagar/casebase/blob/master/references/Saarela-2015-Lifetime_Data_Analysis.pdf).

First example
-------------

The first example we discuss uses the well-known `veteran` dataset,
which is part of the `survival` package. As we can see below, there is
almost no censoring, and therefore we can get a good visual
representation of the survival function:

``` {.r}
library(survival)
data(veteran)
table(veteran$status)
```

    ## 
    ##   0   1 
    ##   9 128

``` {.r}
evtimes <- veteran$time[veteran$status == 1]
hist(evtimes, nclass = 30, main = '', xlab = 'Survival time (days)', 
     col = 'gray90', probability = TRUE)
tgrid <- seq(0, 1000, by = 10)
lines(tgrid, dexp(tgrid, rate = 1.0/mean(evtimes)), 
      lwd = 2, lty = 2, col = 'red')
```

![](smoothHazard_files/figure-markdown/unnamed-chunk-1-1.png)

As we can see, the empirical survival function ressembles an exponential
distribution.

We will first try to estimate the hazard function parametrically using
some well-known regression routines. But first, we will reformat the
data slightly.

``` {.r}
veteran$prior <- factor(veteran$prior, levels = c(0, 10), labels = c("no","yes"))
veteran$celltype <- factor(veteran$celltype, 
                           levels = c('large', 'squamous', 'smallcell', 'adeno'))
veteran$trt <- factor(veteran$trt, levels = c(1, 2), labels = c("standard", "test"))
```

Using the `eha` package, we can fit a Weibull form, with different
values of the shape parameter. For `shape = 1`, we get an exponential
distribution:

``` {.r}
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
    ##               no    0.653     0         1           (reference)
    ##              yes    0.347     0.049     1.051     0.227     0.827 
    ## celltype 
    ##            large    0.269     0         1           (reference)
    ##         squamous    0.421    -0.377     0.686     0.273     0.166 
    ##        smallcell    0.206     0.443     1.557     0.261     0.090 
    ##            adeno    0.104     0.736     2.087     0.294     0.012 
    ## trt 
    ##         standard    0.477     0         1           (reference)
    ##             test    0.523     0.220     1.246     0.199     0.269 
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

If we take `shape = 0`, the shape parameter is estimated along with the
regression coefficients:

``` {.r}
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
    ##               no    0.653     0         1           (reference)
    ##              yes    0.347     0.047     1.048     0.229     0.836 
    ## celltype 
    ##            large    0.269     0         1           (reference)
    ##         squamous    0.421    -0.428     0.651     0.278     0.123 
    ##        smallcell    0.206     0.462     1.587     0.262     0.078 
    ##            adeno    0.104     0.792     2.208     0.300     0.008 
    ## trt 
    ##         standard    0.477     0         1           (reference)
    ##             test    0.523     0.246     1.279     0.203     0.224 
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

``` {.r}
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
    ## prioryes           7.159e-02  1.074e+00  2.323e-01  0.308  0.75794    
    ## celltypesquamous  -4.013e-01  6.695e-01  2.827e-01 -1.420  0.15574    
    ## celltypesmallcell  4.603e-01  1.584e+00  2.662e-01  1.729  0.08383 .  
    ## celltypeadeno      7.948e-01  2.214e+00  3.029e-01  2.624  0.00869 ** 
    ## trttest            2.946e-01  1.343e+00  2.075e-01  1.419  0.15577    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                   exp(coef) exp(-coef) lower .95 upper .95
    ## karno                0.9677     1.0334    0.9573    0.9782
    ## diagtime             1.0001     0.9999    0.9823    1.0182
    ## age                  0.9913     1.0087    0.9734    1.0096
    ## prioryes             1.0742     0.9309    0.6813    1.6937
    ## celltypesquamous     0.6695     1.4938    0.3847    1.1651
    ## celltypesmallcell    1.5845     0.6311    0.9403    2.6699
    ## celltypeadeno        2.2139     0.4517    1.2228    4.0084
    ## trttest              1.3426     0.7448    0.8939    2.0166
    ## 
    ## Concordance= 0.736  (se = 0.03 )
    ## Rsquare= 0.364   (max possible= 0.999 )
    ## Likelihood ratio test= 62.1  on 8 df,   p=1.799e-10
    ## Wald test            = 62.37  on 8 df,   p=1.596e-10
    ## Score (logrank) test = 66.74  on 8 df,   p=2.186e-11

As we can see, all three models are significant, and they give similar
information: `karno` and `celltype` are significant predictors, both
treatment is not.

The method available in this package makes use of *case-base sampling*.
That is, person-moments are randomly sampled across the entire follow-up
time, with some moments corresponding to cases and others to controls.
By sampling person-moments instead of individuals, we can then use
logistic regression to fit smooth-in-time parametric hazard functions.
See the previous section for more details.

First, we will look at the follow-up time by using population-time
plots:

``` {.r}
# create popTime object
pt_veteran <- casebase::popTime(data = veteran)
```

    ## 'time' will be used as the time variable

    ## 'status' will be used as the event variable

    ## Sampling from all remaining individuals under study,
    ##                     regardless of event status

``` {.r}
class(pt_veteran)
```

    ## [1] "popTime"    "data.table" "data.frame"

``` {.r}
# plot method for objects of class 'popTime'
plot(pt_veteran)
```

![](smoothHazard_files/figure-markdown/unnamed-chunk-6-1.png)

Population-time plots are a useful way of visualizing the total
follow-up experience, where individuals appear on the y-axis, and
follow-up time on the x-axis; each individual's follow-up time is
represented by a gray line segment. For convenience, we have ordered the
patients according to their time-to-event, and each event is represented
by a red dot. The censored observations (of which there is only a few)
correspond to the grey lines which do not end with a red dot.

Next, we use case-base sampling to fit a parametric hazard function via
logistic regression. First, we will include time as a linear term; as
noted above, this corresponds to an Gompertz hazard.

``` {.r}
library(casebase)
model4 <- fitSmoothHazard(status ~ time + karno + diagtime + age + prior +
             celltype + trt, data = veteran, ratio = 100, type = "uniform")
```

    ## 'time' will be used as the time variable

``` {.r}
summary(model4)
```

    ## 
    ## Call:
    ## glm(formula = formula, family = binomial, data = sampleData)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -0.4235  -0.1519  -0.1217  -0.0997   3.3952  
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -3.0248787  0.7163987  -4.222 2.42e-05 ***
    ## time               0.0004681  0.0006438   0.727   0.4672    
    ## karno             -0.0308669  0.0052061  -5.929 3.05e-09 ***
    ## diagtime           0.0029721  0.0092760   0.320   0.7487    
    ## age               -0.0033575  0.0092342  -0.364   0.7162    
    ## prioryes           0.0079221  0.2318095   0.034   0.9727    
    ## celltypesquamous  -0.3966056  0.2831742  -1.401   0.1613    
    ## celltypesmallcell  0.4589966  0.2614418   1.756   0.0792 .  
    ## celltypeadeno      0.7336958  0.3005094   2.442   0.0146 *  
    ## trttest            0.2078887  0.2008673   1.035   0.3007    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1436.2  on 12927  degrees of freedom
    ## Residual deviance: 1367.9  on 12918  degrees of freedom
    ## AIC: 1387.9
    ## 
    ## Number of Fisher Scoring iterations: 8

Since the output object from `fitSmoothHazard` inherits from the `glm`
class, we see a familiar result when using the function `summary`.

The main purpose of fitting smooth hazard functions is that it is then
relatively easy to compute absolute risks. For example, we can use the
function `absoluteRisk` to compute the mean absolute risk at 90 days,
which can then be compared to the empirical measure.

``` {.r}
absoluteRisk(object = model4, time = 90)
```

    ## [1] 0.5724214

``` {.r}
ftime <- veteran$time
mean(ftime <= 90)
```

    ## [1] 0.5547445

We can also fit a Weibull hazard by using a logarithmic term for time:

``` {.r}
model5 <- fitSmoothHazard(status ~ log(time) + karno + diagtime + age + prior +
             celltype + trt, data = veteran, ratio = 100, type = "uniform")
```

    ## 'time' will be used as the time variable

``` {.r}
summary(model5)
```

    ## 
    ## Call:
    ## glm(formula = formula, family = binomial, data = sampleData)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -0.4722  -0.1516  -0.1192  -0.0965   3.4179  
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -3.0476673  0.7437538  -4.098 4.17e-05 ***
    ## log(time)          0.0880720  0.0722409   1.219   0.2228    
    ## karno             -0.0335658  0.0055127  -6.089 1.14e-09 ***
    ## diagtime           0.0009047  0.0092804   0.097   0.9223    
    ## age               -0.0056310  0.0092743  -0.607   0.5437    
    ## prioryes           0.0873113  0.2316784   0.377   0.7063    
    ## celltypesquamous  -0.4263186  0.2785697  -1.530   0.1259    
    ## celltypesmallcell  0.4588368  0.2638831   1.739   0.0821 .  
    ## celltypeadeno      0.8878296  0.3031555   2.929   0.0034 ** 
    ## trttest            0.2187086  0.2031419   1.077   0.2816    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1436.2  on 12927  degrees of freedom
    ## Residual deviance: 1362.1  on 12918  degrees of freedom
    ## AIC: 1382.1
    ## 
    ## Number of Fisher Scoring iterations: 8

With case-base sampling, it is straightforward to fit a semi-parametric
hazard function using splines, which can then be used to estimate the
mean absolute risk.

``` {.r}
# Fit a spline for time
library(splines)
model6 <- fitSmoothHazard(status ~ bs(time) + karno + diagtime + age + prior +
             celltype + trt, data = veteran, ratio = 100, type = "uniform")
```

    ## 'time' will be used as the time variable

``` {.r}
summary(model6)
```

    ## 
    ## Call:
    ## glm(formula = formula, family = binomial, data = sampleData)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -0.4246  -0.1535  -0.1201  -0.0955   3.5224  
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -3.065785   0.726464  -4.220 2.44e-05 ***
    ## bs(time)1          1.660880   1.019394   1.629   0.1033    
    ## bs(time)2         -2.849851   1.768409  -1.612   0.1071    
    ## bs(time)3          1.915546   1.020531   1.877   0.0605 .  
    ## karno             -0.030718   0.005390  -5.699 1.20e-08 ***
    ## diagtime           0.001875   0.009434   0.199   0.8425    
    ## age               -0.006177   0.009318  -0.663   0.5074    
    ## prioryes          -0.017958   0.236079  -0.076   0.9394    
    ## celltypesquamous  -0.395571   0.283274  -1.396   0.1626    
    ## celltypesmallcell  0.486887   0.263709   1.846   0.0648 .  
    ## celltypeadeno      0.752699   0.302629   2.487   0.0129 *  
    ## trttest            0.290720   0.207293   1.402   0.1608    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1436.2  on 12927  degrees of freedom
    ## Residual deviance: 1366.5  on 12916  degrees of freedom
    ## AIC: 1390.5
    ## 
    ## Number of Fisher Scoring iterations: 8

``` {.r}
absoluteRisk(object = model6, time = 90)
```

    ## [1] 0.5740258

As we can see from the summary, there is little evidence that splines
actually improve the fit. Moreover, we can see that estimated individual
absolute risks are essentially the same when using either a linear term
or splines:

``` {.r}
linearRisk <- absoluteRisk(object = model4, time = 90, newdata = veteran)
splineRisk <- absoluteRisk(object = model6, time = 90, newdata = veteran)

plot(linearRisk, splineRisk,
     xlab = "Linear", ylab = "Splines", pch = 19)
abline(a = 0, b = 1, lty = 2, lwd = 2, col = 'red')
```

![](smoothHazard_files/figure-markdown/unnamed-chunk-11-1.png)

These last three models give similar information as the first three,
i.e. the main predictors for the hazard are `karno` and `celltype`, with
treatment being non-significant. Moreover, by explicitely including the
time variable in the formula, we see that it is not significant; this is
evidence that the true hazard is exponential.

Finally, we can look at the estimates of the coefficients for the Cox
model, as well as the last three models (CB stands for "case-base"):

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
-0.0309
</td>
<td style="text-align:right;">
-0.0336
</td>
<td style="text-align:right;">
-0.0307
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
0.0030
</td>
<td style="text-align:right;">
0.0009
</td>
<td style="text-align:right;">
0.0019
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
-0.0034
</td>
<td style="text-align:right;">
-0.0056
</td>
<td style="text-align:right;">
-0.0062
</td>
</tr>
<tr>
<td style="text-align:left;">
prioryes
</td>
<td style="text-align:right;">
0.0716
</td>
<td style="text-align:right;">
0.0079
</td>
<td style="text-align:right;">
0.0873
</td>
<td style="text-align:right;">
-0.0180
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
-0.3966
</td>
<td style="text-align:right;">
-0.4263
</td>
<td style="text-align:right;">
-0.3956
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
0.4590
</td>
<td style="text-align:right;">
0.4588
</td>
<td style="text-align:right;">
0.4869
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
0.7337
</td>
<td style="text-align:right;">
0.8878
</td>
<td style="text-align:right;">
0.7527
</td>
</tr>
<tr>
<td style="text-align:left;">
trttest
</td>
<td style="text-align:right;">
0.2946
</td>
<td style="text-align:right;">
0.2079
</td>
<td style="text-align:right;">
0.2187
</td>
<td style="text-align:right;">
0.2907
</td>
</tr>
</tbody>
</table>
Cumulative Incidence Curves
---------------------------

Here we show how to calculate the cumulative incidence curves for a
specific risk profile using the following equation:

$$ CI(x, t) = 1 - exp\left[ - \int_0^t h(x, u) \textrm{d}u \right] $$
where \\( h(x, t) \\) is the hazard function, \\( t \\) denotes the
numerical value (number of units) of a point in prognostic/prospective
time and \\( x \\) is the realization of the vector \\( X \\) of
variates based on the patient's profile and intervention (if any).

We compare the cumulative incidence functions from the fully-parametric
fit using case base sampling, with those from the Cox model:

``` {.r}
# define a specific covariate profile
new_data <- data.frame(trt = "test", 
                       celltype = "adeno", 
                       karno = median(veteran$karno), 
                       diagtime = median(veteran$diagtime),
                       age = median(veteran$age),
                       prior = "no")

# calculate cumulative incidence using casebase model
smooth_risk <- absoluteRisk(object = model4, time = seq(0,300, 1), 
                            newdata = new_data)

# cumulative incidence function for the Cox model
plot(survfit(model3, newdata=new_data),
     xlab = "Days", ylab="Cumulative Incidence (%)", fun = "event",
     xlim = c(0,300), conf.int = F, col = "red", 
     main = sprintf("Estimated Cumulative Incidence (risk) of Lung Cancer\ntrt = test, celltype = adeno, karno = %g,\ndiagtime = %g, age = %g, prior = no", median(veteran$karno), median(veteran$diagtime), 
                    median(veteran$age)))

# add casebase curve with legend
lines(seq(0,300, 1), smooth_risk[1,], type = "l", col = "blue")
legend("bottomright", 
       legend = c("semi-parametric (Cox)", "parametric (casebase)"), 
       col = c("red","blue"),
       lty = c(1, 1), 
       bg = "gray90")
```

![](smoothHazard_files/figure-markdown/unnamed-chunk-13-1.png)

Session information
-------------------

    ## R version 3.3.2 (2016-10-31)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.1 LTS
    ## 
    ## attached base packages:
    ## [1] splines   stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] casebase_0.1.0  eha_2.4-4       survival_2.40-1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.10       knitr_1.15.1       magrittr_1.5      
    ##  [4] munsell_0.4.3      colorspace_1.3-2   lattice_0.20-34   
    ##  [7] highr_0.6          stringr_1.2.0      plyr_1.8.4        
    ## [10] tools_3.3.2        grid_3.3.2         data.table_1.10.4 
    ## [13] gtable_0.2.0       pacman_0.4.1       htmltools_0.3.6   
    ## [16] assertthat_0.1     lazyeval_0.2.0     yaml_2.1.14       
    ## [19] rprojroot_1.2      digest_0.6.12      tibble_1.2        
    ## [22] Matrix_1.2-7.1     ggplot2_2.2.1      VGAM_1.0-3        
    ## [25] evaluate_0.10      rmarkdown_1.3.9003 labeling_0.3      
    ## [28] stringi_1.1.3      scales_0.4.1       backports_1.0.5   
    ## [31] stats4_3.3.2
