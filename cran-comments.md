## Reason for submission

CRAN check results for casebase 0.10.6 showed ERRORs on most platforms, caused
by a breaking change in visreg 3.0 (print.cond renamed to print_cond, and
main/xlab/ylab/gg no longer accepted via ...). This release fixes
plot.singleEventCB() and the affected vignette/tests accordingly. It also
fixes a NOTE for a dead URL in man/support.Rd (the referenced Vanderbilt
Biostatistics wiki page has been retired).

## R CMD check results

There were no ERRORs or WARNINGs

## Test environments

* ubuntu-latest (on GH-Actions), R-devel, R-release, R-oldrel-1
* Windows Server (on GH-Actions), R-release
* macOS (on GH-Actions), R-release
* rocker/verse Docker image (R 4.6.1) with visreg 3.0.0, local
* Win-builder, R-devel and R-release
* mac-builder, R-release

## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTES.

## Reverse dependencies

casebase has one reverse dependency, riskRegression, which only calls
casebase::absoluteRisk(). This release does not modify absoluteRisk() or any
of its dependencies -- the changes are scoped to plot.singleEventCB() (and
its visreg integration) and package documentation, so riskRegression is
unaffected by this release.
