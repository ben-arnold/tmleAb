# tmleAb
tmleAb: Targeted maximum likelihood estimation for antibody data


## Description
Targeted maximum likelihood estimation for antibody measurements in R. In addition to estimating age-dependent antibody curves, the `tmleAb` package makes it easy to estimate marginally adjusted means from the curves, and compare marginally adjusted means between populations. This package is first and foremost a convenience wrapper: the heavy lifting for machine learning is conducted using the `SuperLearner` package, and statistical estimation is conducted using the `tmle` package. However, it also provides some methods for parametric models applied to antibody data (such as antibody acquisition models). 

## Installation

You can install this package easily using the `devtools` package in R. 
```
library(devtools)
install_github("ben-arnold/tmleAb")
```

At this point, the package is fully functional but should be considered a beta version. Check back here for updates: we have many updated features planned and will eventually distribute the package through CRAN.  If you use the package and have any comments please send them along -- we'd love to hear from you (email: benarnold@berkeley.edu).
