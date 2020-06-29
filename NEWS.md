# casebase 0.9.0

* *Major new release*
* Some minor bug fixes
* Trapezoidal rule to perform numerical integration for absolute risk estimation, providing significant speed up
* Major refactoring of absoluteRisk
* Users now have further control on the output of `absoluteRisk` using the arguments `type` and `addZero`.
* The output of `absoluteRisk` now always contains the time variable in the first column, regardless of the length of `time`. This will break earlier code that depended on the previous behaviour.
* Population time plots now use `geom_ribbon` instead of `geom_segment`
* Population time functions now allowed for more flexible plots with user defined arguments including sequentially adding base, case, and competing event series
* New plotting method for time-dependent hazard functions and hazard ratios. These include confidence intervals. See `plot.singleEventCB`. The hazard function plot requires `visreg` package to be installed. 
* New plotting method for cumulative incidence and survival curves. See `plot.absRiskCB`.
* No longer importing the entire namespace of `data.table` and `ggplot2`. 

# casebase 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* First release of the casebase package


