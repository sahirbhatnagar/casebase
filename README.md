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


## Citation

To cite `casebase` in publications, please use

```R
citation('casebase')
```

```
Bhatnagar S, Turgeon M, Saarela O and Hanley J (2017). 
casebase: Fitting Flexible Smooth-in-Time Hazards and Risk Functions via Logistic and Multinomial Regression. 
R package version 0.1.0, <URL:https://CRAN.R-project.org/package=casebase>.

Hanley, James A., and Olli S. Miettinen. 
Fitting smooth-in-time prognostic risk functions via logistic regression. 
International Journal of Biostatistics 5.1 (2009): 1125-1125.

Saarela, Olli. A case-base sampling method for estimating recurrent event intensities. 
Lifetime data analysis 22.4 (2016): 589-605.

If competing risks analyis is used, please also cite:

Saarela, Olli, and Elja Arjas. Non-parametric Bayesian Hazard Regression for Chronic Disease Risk Assessment. 
Scandinavian Journal of Statistics 42.2 (2015): 609-626.
```

For BibTeX users:

```R
toBibtex(citation('casebase'))
```

```
@Manual{casebase-package,
  title = {casebase: Fitting Flexible Smooth-in-Time Hazards and Risk Functions via Logistic and Multinomial Regression},
  author = {Sahir Bhatnagar and Maxime Turgeon and Olli Saarela and James Hanley},
  year = {2017},
  note = {R package version 0.1.0},
  url = {https://CRAN.R-project.org/package=casebase},
}

@Article{,
  title = {Fitting smooth-in-time prognostic risk functions via logistic regression},
  author = {James A Hanley and Olli S Miettinen},
  journal = {International Journal of Biostatistics},
  volume = {5},
  number = {1},
  pages = {1125--1125},
  year = {2009},
  publisher = {Berkeley Electronic Press},
}

@Article{,
  title = {A case-base sampling method for estimating recurrent event intensities},
  author = {Olli Saarela},
  journal = {Lifetime data analysis},
  volume = {22},
  number = {4},
  pages = {589--605},
  year = {2016},
  publisher = {Springer},
}

@Article{,
  title = {Non-parametric Bayesian Hazard Regression for Chronic Disease Risk Assessment},
  author = {Olli Saarela and Elja Arjas},
  journal = {Scandinavian Journal of Statistics},
  year = {2015},
  volume = {42},
  number = {2},
  pages = {609--626},
  publisher = {Wiley Online Library},
}
```

## References

<ol>
<li>
<p>Hanley, James A, and Olli S Miettinen. 2009. "Fitting Smooth-in-Time Prognostic Risk Functions via Logistic Regression." <em>The International Journal of Biostatistics</em> 5 (1).</p>
</li>
<li>
<p>Saarela, Olli, and Elja Arjas. 2015. "Non-Parametric Bayesian Hazard Regression for Chronic Disease Risk Assessment." <em>Scandinavian Journal of Statistics</em> 42 (2). Wiley Online Library: 609–26.</p>
</li>
<li>
<p>Saarela, Olli. 2015. "A Case-Base Sampling Method for Estimating Recurrent Event Intensities." <em>Lifetime Data Analysis</em>. Springer, 1–17.</p>
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
