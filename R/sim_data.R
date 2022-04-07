#' Simulate data with measurment error from a simple linear regression
#'
#' @param n_sim Number of data points to simulate
#' @param min_x Minimum x value
#' @param max_x Maximum x value
#' @param alpha Regression intercept
#' @param beta  Regression slope
#' @param sigma Nugget variation
#' @param x_err x measurement error
#' @param y_err y measurement error
#'
#' @return Simulated dataset with columns x, x_err, y, y_err
#' @export
#'
#' @examples
#' sim_reg(n_sim = 50)
sim_reg <- function(n_sim = 50,
                    min_x = 0,
                    max_x = 2,
                    alpha = 0,
                    beta = 1,
                    sigma = 0.1,
                    x_err = 0.1,
                    y_err = 0.1) {
  sdobs <- x_err
  taux <- 1 / (sdobs * sdobs)
  truex <- stats::runif(n_sim, min_x, max_x)
  obsx <- stats::rnorm(n_sim, truex, sdobs)
  errory <- stats::rnorm(n_sim, 0, y_err)
  true_reg <- alpha + beta * truex
  obsy <- true_reg + errory
  parms <- data.frame(alpha, beta)

  dat <- tibble::tibble(
    x = obsx,
    x_err = x_err,
    y = obsy,
    y_err = y_err,
    true_line = true_reg,
    true_x = truex
  )

  return(dat)
}


#' Simulate data with measurment error from a change-point linear regression
#'
#' @param n_sim Number of data points to simulate
#' @param min_x Minimum x value
#' @param max_x Maximum x value
#' @param alpha Regression intercept
#' @param beta  Regression slope
#' @param cp Change point
#' @param sigma Nugget variation
#' @param x_err x measurement error
#' @param y_err y measurement error
#'
#' @return Simulated dataset with columns x, x_err, y, y_err
#' @export
#'
#' @examples
#' sim_cp(n_sim = 50)
sim_cp <- function(n_sim = 50,
                   min_x = 0,
                   max_x = 2,
                   alpha = 0,
                   beta = c(1, 2),
                   cp = 1,
                   sigma = 0.1,
                   x_err = 0.1,
                   y_err = 0.1) {
  sdobs <- x_err
  taux <- 1 / (sdobs * sdobs)
  truex <- stats::runif(n_sim, min_x, max_x)
  beta_index <- ifelse(truex < 1, 1, 2)
  obsx <- stats::rnorm(n_sim, truex, sdobs)
  errory <- stats::rnorm(n_sim, 0, y_err)
  true_cp <- rep(NA, n_sim)
  for (i in 1:n_sim) true_cp[i] <- alpha + beta[beta_index[i]] * (truex[i] - cp)

  obsy <- true_cp + errory

  dat <- tibble::tibble(
    x = obsx,
    x_err = x_err,
    y = obsy,
    y_err = y_err,
    true_line = true_cp,
    true_x = truex
  )

  return(dat)
}

#' Simulate data with measurment error from a Gaussian process regression
#'
#' @param n_sim number of data points to simulate
#' @param min_x Minimum x value
#' @param max_x Maximum x value
#' @param alpha regression intercept
#' @param sigma_g GP variance parameter
#' @param phi GP correlation parameter
#' @param sigma nugget variation
#' @param x_err x measurement error
#' @param y_err y measurement error
#'
#' @return Simulated dataset with columns x, x_err, y, y_err
#' @export
#'
#' @examples
#' sim_gp(n_sim = 50)
sim_gp <- function(n_sim = 50,
                   min_x = 0,
                   max_x = 2,
                   alpha = 0,
                   sigma_g = 2,
                   phi = 2,
                   sigma = 0.1,
                   x_err = 0.1,
                   y_err = 0.1) {
  true_x <- sort(stats::runif(min_x, max_x, n = n_sim))
  n_obs <- n_sim
  x_star <- true_x

  ### Predictive distribution for the GP
  Sigma <- (sigma_g^2) * exp(-(phi^2) * fields::rdist(true_x)^2)
  gp <- c(mvtnorm::rmvnorm(n = 1, sigma = Sigma))


  sim_eiv <- "
data{
   for(i in 1:n_sim)
{
    y[i] ~ dnorm(alpha + gp[i],(sigma + y_err)^-2)
    x[i] ~ dnorm(mu_x[i],x_err^-2)
  }
  }
model{
fake <- 0
}
"

  data <- list(
    n_sim = n_sim,
    mu_x = true_x,
    alpha = alpha,
    sigma = sigma,
    x_err = x_err,
    y_err = y_err,
    gp = gp
  )

  out <- runjags::run.jags(sim_eiv,
    data = data,
    monitor = c("y", "x"),
    sample = 1,
    n.chains = 1,
    summarise = FALSE
  )


  sim_dat <- coda::as.mcmc(out)

  x <- as.vector(sim_dat)[(n_sim + 1):(n_sim * 2)]
  y <- as.vector(sim_dat)[1:n_sim]

  dat <- tibble::tibble(
    x = x,
    x_err = x_err,
    y = y,
    y_err = y_err,
    true_line = gp,
    true_x = true_x
  )

  return(dat)
}
