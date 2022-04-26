# EIVmodels 

In statistics, errors-in-variables (EIV) models, or measurement error models, are  models that account for measurement errors in both the independent (predictor) and dependent (outcome) variables. EIVmodels is a R package designed specifically to account for measurement errors within some commonly used models (linear regression, change-point regression, (Integrated) Gaussian process regression) when analysing time-dependent data derived from paleoenvironmental reconstructions.The models are implemented in a Bayesian framework using the JAGS (Just Another Gibbs Sampler) software.

## Requirements

  - The package requires the installation of the JAGS software. Click to [download JAGS](https://sourceforge.net/projects/mcmc-jags/).

## Getting started

See [Vignettes](https://ncahill89.github.io/EIVmodels/articles/eivmodels.html).

## Installation

  - This package is not currently on cran so you can download from Github. Make sure to have the `devtools` package installed and then execute the following: 

```
devtools::install_github("ncahill89/EIVmodels")
```
