
[![Build Status](https://travis-ci.org/ben-arnold/tmleAb.svg)](https://travis-ci.org/ben-arnold/tmleAb)
![CRAN](http://www.r-pkg.org/badges/version/tmleAb)

tmleAb
======

tmleAb: Targeted maximum likelihood estimation for antibody data

The `tmleAb` package estimates age-dependent antibody curves with ensemble machine learning, makes it easier to estimate marginally adjusted means from the curves, and compare marginally adjusted means between populations.

This package is mostly a convenience wrapper: the heavy lifting for machine learning is conducted using the `SuperLearner` package, and statistical estimation is conducted using the `tmle` package. In addition to these features, `tmleAb` also provides some methods for parametric models applied to antibody data (such as antibody acquisition models).

The package's web interface is: <http://ben-arnold.github.io/tmleAb/>

An article that describes many of the methods implmented in the package is here: <http://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005616>

With an associated GitHub repository here: <https://github.com/ben-arnold/antibody-curves>

Installation
------------

Install `tmleAb` using the `devtools` package in R.

    library(devtools)
    install_github("ben-arnold/tmleAb")

The `tmleAb` package is in occasional development, so periodically check for updates. The package is not intended to be commercial-grade software, but it is intended to make it easier to use ensemble machine learning and targeted maximum likelihood estimation for the analysis of antibody measurements.  If you use the package and run into issues and/or have comments, please write. 

The package vignette (see `Articles` on the top of the package's website) includes some details about the package and a few examples of how to use its main features.
