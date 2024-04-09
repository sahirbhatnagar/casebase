# riskRegression

<details>

* Version: 2023.12.21
* GitHub: https://github.com/tagteam/riskRegression
* Source code: https://github.com/cran/riskRegression
* Date/Publication: 2023-12-19 17:00:02 UTC
* Number of recursive dependencies: 186

Run `revdepcheck::revdep_details(, "riskRegression")` for more info

</details>

## In both

*   checking whether package ‘riskRegression’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/maxturgeon/git_repos/casebase/revdep/checks.noindex/riskRegression/new/riskRegression.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘riskRegression’ ...
** package ‘riskRegression’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
using C++ compiler: ‘Homebrew clang version 14.0.6’
using SDK: ‘MacOSX14.4.sdk’
/opt/homebrew/opt/llvm/bin/clang++ -fopenmp -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/maxturgeon/git_repos/casebase/revdep/library.noindex/riskRegression/Rcpp/include' -I'/Users/maxturgeon/git_repos/casebase/revdep/library.noindex/riskRegression/RcppArmadillo/include' -I/opt/homebrew/opt/gettext/include -I/opt/homebrew/opt/llvm/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include    -fPIC  -g -O3 -Wall -pedantic -std=c++11 -mtune=native -pipe -c IC-Nelson-Aalen-cens-time.cpp -o IC-Nelson-Aalen-cens-time.o
In file included from IC-Nelson-Aalen-cens-time.cpp:2:
In file included from ./IC-Nelson-Aalen-cens-time.h:1:
In file included from /Users/maxturgeon/git_repos/casebase/revdep/library.noindex/riskRegression/RcppArmadillo/include/RcppArmadillo.h:29:
...
        return complex<_Tp>(abs(__x.real()), __x.imag());
                            ^
/opt/homebrew/opt/llvm/bin/../include/c++/v1/cmath:338:1: note: using declaration annotated with 'using_if_exists' here
using ::abs _LIBCPP_USING_IF_EXISTS;
^
fatal error: too many errors emitted, stopping now [-ferror-limit=]
20 errors generated.
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
using C++ compiler: ‘Homebrew clang version 14.0.6’
using SDK: ‘MacOSX14.4.sdk’
/opt/homebrew/opt/llvm/bin/clang++ -fopenmp -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/maxturgeon/git_repos/casebase/revdep/library.noindex/riskRegression/Rcpp/include' -I'/Users/maxturgeon/git_repos/casebase/revdep/library.noindex/riskRegression/RcppArmadillo/include' -I/opt/homebrew/opt/gettext/include -I/opt/homebrew/opt/llvm/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include    -fPIC  -g -O3 -Wall -pedantic -std=c++11 -mtune=native -pipe -c IC-Nelson-Aalen-cens-time.cpp -o IC-Nelson-Aalen-cens-time.o
In file included from IC-Nelson-Aalen-cens-time.cpp:2:
In file included from ./IC-Nelson-Aalen-cens-time.h:1:
In file included from /Users/maxturgeon/git_repos/casebase/revdep/library.noindex/riskRegression/RcppArmadillo/include/RcppArmadillo.h:29:
...
        return complex<_Tp>(abs(__x.real()), __x.imag());
                            ^
/opt/homebrew/opt/llvm/bin/../include/c++/v1/cmath:338:1: note: using declaration annotated with 'using_if_exists' here
using ::abs _LIBCPP_USING_IF_EXISTS;
^
fatal error: too many errors emitted, stopping now [-ferror-limit=]
20 errors generated.
make: *** [IC-Nelson-Aalen-cens-time.o] Error 1
ERROR: compilation failed for package ‘riskRegression’
* removing ‘/Users/maxturgeon/git_repos/casebase/revdep/checks.noindex/riskRegression/old/riskRegression.Rcheck/riskRegression’


```
