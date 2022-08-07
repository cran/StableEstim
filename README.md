<!-- badges: start -->
  [![CRANStatusBadge](http://www.r-pkg.org/badges/version/StableEstim)](https://cran.r-project.org/package=StableEstim)
  [![R-CMD-check](https://github.com/GeoBosh/StableEstim/workflows/R-CMD-check/badge.svg)](https://github.com/GeoBosh/StableEstim/actions)
<!-- badges: end -->


# Installing StableEstim

The [latest stable version](https://cran.r-project.org/package=StableEstim) is on
CRAN.

    install.packages("StableEstim")

You can install the [development version](https://github.com/GeoBosh/StableEstim) of
`StableEstim` from Github:

    library(devtools)
    install_github("GeoBosh/StableEstim")


# Overview

A collection of methods to estimate the four parameters of stable
distributions. The package also provides functions to compute
characteristic functions and tools to run Monte Carlo simulations.

The main functions of package `StableEstim` are briefly described below:


* main function: `Estim()` estimates the parameters by various
  methods. Also gives the associated asymptotic properties of the
  estimators.

* estimation functions for specific methods: these functions are called by `Estim()` but can be used directly, as well. The methods provided so far are:

  - the maximum-likelihood (`MLParametersEstim()`),

  - the generalised method of moments with a finite (`GMMParametersEstim()`)
    or continuum moment conditions (`CgmmParametersEstim()`),

  - the iterative Koutrouvelis regression method
    (`KoutParametersEstim()`),

  - the fast Kogon-McCulloch method used for first guess estimation
    (`IGParametersEstim`).
      
* characteristic function: the characteristic function (`ComplexCF()`)
  and its Jacobian (`jacobianComplexCF()`) can be computed and will
  return a vector (respectively a matrix) of complex numbers.

* Monte Carlo simulation: a tool to run a Monte Carlo simulation
  (`Estim_Simulation()`) is provided and can save output files and/or
  produce statistical summary.



The package is developed by Tarak Kharrat and Georgi N.Boshnakov.
