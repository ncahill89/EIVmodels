# EIVmodels 


EIVmodels is a R package designed to run some commonly used models (linear regression, change-point regression, (Integrated) Gaussian process regression) accounting for measurement error in a Bayesian framework. Learn more in `vignette("eivmodels")`.

## Installation

You will need to install the JAGS (Just Another Gibbs Sampler) software in order for the functions in the package to work. You can download JAGS from [here](https://sourceforge.net/projects/mcmc-jags/).

This package is not currently on cran so you can download from Github. Make sure to have the `devtools` package installed and the execute the following: 

```
# Install development version from Github
devtools::install_github("ncahill89/EIVmodels")
```


