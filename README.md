# tmleAb
tmleAb: Targeted maximum likelihood estimation for antibody data


In addition to estimating age-dependent antibody curves, the `tmleAb` package makes it easy to estimate marginally adjusted means from the curves, and compare marginally adjusted means between populations. 

This package is first and foremost a convenience wrapper: the heavy lifting for machine learning is conducted using the `SuperLearner` package, and statistical estimation is conducted using the `tmle` package. In addition to these features, the package also provides some methods for parametric models applied to antibody data (such as antibody acquisition models).

The package's more convenient web interface is: http://ben-arnold.github.io/tmleAb/

## Installation

You can install this package easily using the `devtools` package in R. 
```
library(devtools)
install_github("ben-arnold/tmleAb")
```

The `tmleAb` package is a development version. Check back here for updates: we have many updated features planned.  If you use the package and have any comments please send them along -- we'd love to hear from you.  The package vignette includes some details about the package and a few examples of how to use its key features. 


