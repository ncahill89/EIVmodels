# EIVmodels 

In statistics, errors-in-variables (EIV) models, or measurement error models, are regression models that account for measurement errors in both the independent (predictor) and dependent (outcome) variables. EIVmodels is a R package designed to account for measurement errors in some commonly used models (linear regression, change-point regression, (Integrated) Gaussian process regression). The models are implemented in a Bayesian framework using the JAGS (Just Another Gibbs Sampler) software.

## Requirements

  - The package requires the installation of the JAGS software. You can download JAGS from [here](https://sourceforge.net/projects/mcmc-jags/).

## Getting started

See Slides, Vignettes and Tutorials

## Installation

  - This package is not currently on cran so you can download from Github. Make sure to have the `devtools` package installed and then execute the following: 

```
devtools::install_github("ncahill89/EIVmodels")
```
