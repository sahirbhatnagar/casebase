# casebase 0.10.3

* Fixed broken links in vignette.

# casebase 0.10.2

* Fixed noLD checks issue as reported by CRAN (Issue 156).

# casebase 0.10.1

* Fixed issue 143 and return the data invisibly with `plot.singleEventCB()` when `type = "hr"`.
* Removed `family = "gbm"` as it wasn't properly tested.
* Added `confint.singleEventCB` to compute confidence bands for the risk (or survival) function.
* Updated `ERSPC` data so that the exposure variable is categorical. This may break previous code explicitly making this conversion, or somehow relying on the numerical coding.

# casebase 0.9.1

* Fixed issue with `plot.singleEventCB()` when `visreg` package is not loaded.
* Improved error message when using `family = "glmnet"` with a single covariate.
* Introduced `summary` method for objects of class `singleEventCB`, and improved the output of `print` by displaying the appropriate function call.

# casebase 0.9.0

This is a *Major new release*

## Breaking changes

* The output of `absoluteRisk()` now always contains the time variable in the first column, regardless of the length of `time`. This will break earlier code that depended on the previous behaviour.
* Population time plots now use `ggplot2::geom_ribbon()` instead of `ggplot2::geom_segment()`. 
* Population time functions now allow for more flexible plots with user defined arguments including sequentially adding base, case, and competing event series. These are now passed as a list to the `*.params` arguments. Several arguments are now deprecated. 
* Removed `popTimeExposure` class and the corresponding `plot` method. `popTime()` now returns an `exposure` attribute which contains the name of the exposure variable in the dataset. The plot method for objects of class `popTime` will use this exposure attribute to create exposure stratified population time plots.

## New features

* Major refactoring of `absoluteRisk()`. Trapezoidal rule to perform numerical integration for absolute risk estimation, providing significant speed up.
* Users now have further control on the output of `absoluteRisk()` using the arguments `type` and `addZero`.
* New plotting method for time-dependent hazard functions and hazard ratios. These include confidence intervals. See `plot.singleEventCB()`. The hazard function plot requires the `visreg` package to be installed. 
* New plotting method for cumulative incidence and survival curves. See `plot.absRiskCB()`.
* When `time` is unspecified, `absoluteRisk()` now computes the cumulative incidence at `ntimes` equidistant points between 0 and the max failure time.
* `absoluteRisk()` can now compute the cumulative incidence for a `"typical"` covariate profile with `newdata = "typical"`. "Typical" corresponds to the median for continuous variables and the mode for factors (each variable is summarised independently).
* Added `eprchd`, `brcancer`, `support` and `simdat` datasets to the package. 
* Implemented `riskRegression::predictRisk()` method for `singleEventCB` objects.
 
## Minor improvements and fixes

* No longer importing the entire namespace of `data.table` and `ggplot2`. 
* Moved from make the docs to pkgdown for package website.
* A warning is given when `family="gbm"` and nonlinear functions of time or interactions are specified.
* Add `singleEventCB` class to object returned by `fitSmoothHazard()`
* Add `absRiskCB` class to object returned by `absoluteRisk()`
* Use `glmnet::prepareX` to convert factors into indicator variables


# casebase 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* First release of the casebase package


