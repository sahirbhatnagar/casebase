# riskRegression

<details>

* Version: 2021.10.10
* GitHub: https://github.com/tagteam/riskRegression
* Source code: https://github.com/cran/riskRegression
* Date/Publication: 2021-10-11 10:30:02 UTC
* Number of recursive dependencies: 175

Run `revdep_details(, "riskRegression")` for more info

</details>

## In both

*   checking whether package ‘riskRegression’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/maxturgeon/git_repos/casebase/revdep/checks.noindex/riskRegression/new/riskRegression.Rcheck/00install.out’ for details.
    ```

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'smcfcs', 'randomForestSRC', 'scam', 'SuperLearner'
    ```

## Installation

### Devel

```
* installing *source* package ‘riskRegression’ ...
** package ‘riskRegression’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.1/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.1/Resources/library/RcppArmadillo/include' -I/usr/local/opt/gettext/include -I/usr/local/opt/llvm/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include   -fPIC  -Wall -g -O2  -c AUCijFun.cpp -o AUCijFun.o
clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.1/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.1/Resources/library/RcppArmadillo/include' -I/usr/local/opt/gettext/include -I/usr/local/opt/llvm/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include   -fPIC  -Wall -g -O2  -c IC-Nelson-Aalen-cens-time.cpp -o IC-Nelson-Aalen-cens-time.o
error: unable to rename temporary 'IC-Nelson-Aalen-cens-time-de0b3997.o.tmp' to output file 'IC-Nelson-Aalen-cens-time.o': 'No space left on device'
1 error generated.
make: *** [IC-Nelson-Aalen-cens-time.o] Error 1
ERROR: compilation failed for package ‘riskRegression’
* removing ‘/Users/maxturgeon/git_repos/casebase/revdep/checks.noindex/riskRegression/new/riskRegression.Rcheck/riskRegression’


```
### CRAN

```
* installing *source* package ‘riskRegression’ ...
** package ‘riskRegression’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.1/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.1/Resources/library/RcppArmadillo/include' -I/usr/local/opt/gettext/include -I/usr/local/opt/llvm/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include   -fPIC  -Wall -g -O2  -c AUCijFun.cpp -o AUCijFun.o
clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.1/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.1/Resources/library/RcppArmadillo/include' -I/usr/local/opt/gettext/include -I/usr/local/opt/llvm/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include   -fPIC  -Wall -g -O2  -c IC-Nelson-Aalen-cens-time.cpp -o IC-Nelson-Aalen-cens-time.o
fatal error: error in backend: IO failure on output stream: No space left on device
make: *** [IC-Nelson-Aalen-cens-time.o] Error 1
ERROR: compilation failed for package ‘riskRegression’
* removing ‘/Users/maxturgeon/git_repos/casebase/revdep/checks.noindex/riskRegression/old/riskRegression.Rcheck/riskRegression’


```
