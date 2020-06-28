# casebase 1.0.0

* *Major new release*
* Some minor bug fixes
* Trapezoidal rule to perform numerical integration for absolute risk estimation, providing significant speed up
* Major refactoring of absoluteRisk
* Users now have further control on the output of `absoluteRisk` using the arguments `type` and `addZero`.
* The output of `absoluteRisk` now always contains the time variable in the first column, regardless of the length of `time`. This will break earlier code that depended on the previous behaviour.
* Population time plots now use `geom_ribbon` instead of `geom_segment`
* Population time functions now allowed for more flexible plots with user defined arguments

# casebase 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* First release of the casebase package


