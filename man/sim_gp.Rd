% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_data.R
\name{sim_gp}
\alias{sim_gp}
\title{Simulate data with measurment error from a Gaussian process regression}
\usage{
sim_gp(
  n_sim = 50,
  min_x = 0,
  max_x = 2,
  alpha = 0,
  sigma_g = 2,
  phi = 2,
  sigma = 0.1,
  x_err = 0.1,
  y_err = 0.1
)
}
\arguments{
\item{n_sim}{number of data points to simulate}

\item{min_x}{Minimum x value}

\item{max_x}{Maximum x value}

\item{alpha}{regression intercept}

\item{sigma_g}{GP variance parameter}

\item{phi}{GP correlation parameter}

\item{sigma}{nugget variation}

\item{x_err}{x measurement error}

\item{y_err}{y measurement error}
}
\value{
Simulated dataset with columns x, x_err, y, y_err
}
\description{
Simulate data with measurment error from a Gaussian process regression
}
\examples{
sim_gp(n_sim = 50)
}
