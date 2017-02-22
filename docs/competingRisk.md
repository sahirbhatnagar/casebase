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

    ##   Length    Class     Mode 
    ##        1 CompRisk       S4

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

    ##   Length    Class     Mode 
    ##        1 CompRisk       S4

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

    ##   Length    Class     Mode 
    ##        1 CompRisk       S4

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

    ## R version 3.3.2 (2016-10-31)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.1 LTS
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] casebase_0.1.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.9        knitr_1.15.1       magrittr_1.5      
    ##  [4] splines_3.3.2      munsell_0.4.3      lattice_0.20-34   
    ##  [7] colorspace_1.3-2   stringr_1.1.0      plyr_1.8.4        
    ## [10] tools_3.3.2        grid_3.3.2         data.table_1.10.4 
    ## [13] gtable_0.2.0       htmltools_0.3.6    survival_2.40-1   
    ## [16] yaml_2.1.14        lazyeval_0.2.0     rprojroot_1.2     
    ## [19] digest_0.6.12      assertthat_0.1     tibble_1.2        
    ## [22] Matrix_1.2-7.1     ggplot2_2.2.1      VGAM_1.0-3        
    ## [25] evaluate_0.10      rmarkdown_1.3.9003 stringi_1.1.2     
    ## [28] scales_0.4.1       backports_1.0.5    stats4_3.3.2

References
----------
