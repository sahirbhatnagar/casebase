# casebase

[![Build Status](https://travis-ci.org/sahirbhatnagar/casebase.svg?branch=master)](https://travis-ci.org/sahirbhatnagar/casebase) [![Coverage Status](https://img.shields.io/codecov/c/github/sahirbhatnagar/casebase/master.svg)](https://codecov.io/github/sahirbhatnagar/casebase?branch=master)

An R package for smooth-in-time fitting of parametric hazard functions

## Installation

You can install the development version of `casebase` from [GitHub](https://github.com/sahirbhatnagar/casebase) with:

```R
install.packages("pacman")
pacman::p_install_gh("sahirbhatnagar/casebase")
```

## Vignette

See the [package website](http://sahirbhatnagar.com/casebase/) for example usage of the functions. This includes

1. [Fitting Smooth Hazard Functions](http://sahirbhatnagar.com/casebase/smoothHazard/)
2. [Competing Risks Analysis](http://sahirbhatnagar.com/casebase/competingRisk/)
3. [Population Time Plots](http://sahirbhatnagar.com/casebase/popTime/)

## Credit

This package is makes use of several existing packages including:

* [`VGAM`](https://cran.r-project.org/package=VGAM) for fitting multinomial logistic regression models
* [`survival`](https://cran.r-project.org/package=survival) for survival models
* [`ggplot2`](https://cran.r-project.org/package=ggplot2) for plotting the population time plots


## References


<ol>
<li>
<p>Efron, Bradley. 1977. "The Efficiency of Cox's Likelihood Function for Censored Data." <em>Journal of the American Statistical Association</em> 72 (359). Taylor &amp; Francis Group: 557–65.</p>
</li>
<li>
<p>Hanley, James A, and Olli S Miettinen. 2009. "Fitting Smooth-in-Time Prognostic Risk Functions via Logistic Regression." <em>The International Journal of Biostatistics</em> 5 (1).</p>
</li>
<li>
<p>Mantel, Nathan. 1973. "Synthetic Retrospective Studies and Related Topics." <em>Biometrics</em>. JSTOR, 479–86.</p>
</li>
<li>
<p>Saarela, Olli. 2015. "A Case-Base Sampling Method for Estimating Recurrent Event Intensities." <em>Lifetime Data Analysis</em>. Springer, 1–17.</p>
</li>
<li>
<p>Saarela, Olli, and Elja Arjas. 2015. "Non-Parametric Bayesian Hazard Regression for Chronic Disease Risk Assessment." <em>Scandinavian Journal of Statistics</em> 42 (2). Wiley Online Library: 609–26.</p>
</li>
<li>
<p>Scrucca, L, A Santucci, and F Aversa. 2010. "Regression Modeling of Competing Risk Using R: An in Depth Guide for Clinicians." <em>Bone Marrow Transplantation</em> 45 (9). Nature Publishing Group: 1388–95.</p>
</li>
</ol>


## Contact

* Issues: <https://github.com/sahirbhatnagar/casebase/issues>
* Pull Requests: <https://github.com/sahirbhatnagar/casebase/>
* e-mail: <sahir.bhatnagar@gmail.com>, <maxime.turgeon@mail.mcgill.ca>


## Latest news

You can see the most recent changes to the package in the [NEWS.md file](https://github.com/sahirbhatnagar/casebase/blob/master/NEWS.md)

## Code of Conduct
 
Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
