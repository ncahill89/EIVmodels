#' Fit errors-in-variables models with JAGS
#'
#' @param dat Input data with columns x,x_err,y,y_err
#' @param model The model to run. Choose from slr, cp, gp, igp
#' @param EIV Indicate whether to use EIV framework. Defaults to TRUE
#' @param n_cp Number of change points if model = "cp" is chosen. Can choose from 1,2,3,4.
#' @param igp_smooth Informs prior for the smoothness (correlation) parameter if model = "igp" is chosen. Choose a value between 0 and 1. Closer to 1 will increase smoothness.
#' @param iter MCMC iterations
#' @param burnin MCMC burnin
#' @param thin MCMC thinning
#' @param scale_factor value to divide the predictor (x) by to change the scale. 1 will result in n0 change. 1000 is recommended if x is years.
#'
#' @return List with data, JAGS input data, model output and the name of the model file.
#' @export
#'
#' @examples
#' dat <- sim_slr(n_sim = 30)
#' mod <- run_mod(dat, model = "slr")
run_mod <- function(dat,
                    model = "slr",
                    EIV = TRUE,
                    n_cp = 1,
                    igp_smooth = 0.2,
                    iter = 5000,
                    burnin = 1000,
                    thin = 4,
                    scale_factor = 1) {

  # Simple Linear Regression Model ------------------------------------------
  x <- x_err <- NULL

  dat <- dat %>% dplyr::mutate(
    x_st = x / scale_factor,
    x_err_st = x_err / scale_factor
  )

  myinitial <- NULL

  if (model == "slr" & !EIV) {
    run_model <-
      "model{
  ## Data model loop
  for(j in 1:n_obs)
  {
  y[j]~dnorm(mu_y[j],tau[j])
  mu_y[j] <- alpha + beta*(x[j])
  tau[j] <- (y_err[j]^2 + sigma^2)^-1

  } # end j loop

 ## Priors
 alpha ~ dnorm(0.0,0.01)
 beta ~ dnorm(0.0,0.01)
 sigma ~ dt(0,4^-2,1)T(0,)

 for(i in 1:n_pred) {mu_pred[i] = alpha + beta*x_pred_st[i]}

}# end
"
  }

  if (model == "slr" && EIV) {
    run_model <-
      "model{
  ## Data model loop
  for(j in 1:n_obs)
  {
  y[j] ~ dnorm(mu_y[j],tau[j])
  x[j] ~ dnorm(mu_x[j],x_err[j]^-2)
  mu_x[j] ~ dnorm(0,0.5^-2)
  mu_y[j] <- alpha + beta*(mu_x[j])
  tau[j] <- (y_err[j]^2 + sigma^2)^-1

  } # end j loop

 ## Priors
 alpha ~ dnorm(0.0,0.01)
 beta ~ dnorm(0.0,0.01)
 sigma ~ dt(0,4^-2,1)T(0,)

 for(i in 1:n_pred) {mu_pred[i] = alpha + beta*x_pred_st[i]}

}# end
"
  }

  if (model == "cp" && EIV == FALSE) {
    if (n_cp == 1) {
      run_model <- "

model{
####CP Regression model

  ###Data Loop
    for(j in 1:n_obs)
  {
  y[j]~dnorm(mu_y[j],tau[j])
  C[j] <- 1+step(x[j]-cp)
  mu_y[j] <- alpha + beta[C[j]]*(x[j]-cp)
  tau[j] <- (y_err[j]^2 + sigma^2)^-1
  }

  ##Priors
  alpha[1] ~ dnorm(0.0,0.01)

  beta[1]~dnorm(0.0,0.01)
  beta[2]~dnorm(0.0,0.01)

  sigma ~ dt(0,4^-2,1)T(0,)

  cp ~ dunif(x_min,x_max)

  for(i in 1:n_pred)
  {
  mu_pred[i] <-alpha + beta[Cstar[i]]*(x_pred_st[i]-cp)
  Cstar[i] <- 1+step(x_pred_st[i]-cp)
  }


}"
    }

    if (n_cp == 2) {
      run_model <- "

model{
####CP Regression model

  ###Data Loop
    for(j in 1:n_obs)
  {
  y[j]~dnorm(mu_y[j],tau[j])
  A[j] <- step(x[j]-cp[1])
  C[j] <- 1+step(x[j]-cp[2])
  mu_y[j] <- alpha[C[j]] + beta[A[j] + C[j]]*(x[j]-cp[C[j]])
  tau[j] <- (y_err[j]^2 + sigma^2)^-1
  }

  ##Priors
  alpha[1] ~ dnorm(0.0,0.01)
  alpha[2] ~ dnorm(0.0,0.01)

  beta[1]~dnorm(0.0,0.01)
  beta[2] <- (alpha[2] - alpha[1])/(cp[2]-cp[1])
  beta[3]~dnorm(0.0,0.01)

  sigma ~ dt(0,4^-2,1)T(0,)

  cp.temp[1] ~ dunif(x_min,x_max)
  cp.temp[2] ~ dunif(x_min,x_max)

  cp[1:2]<-sort(cp.temp)

  for(i in 1:n_pred)
  {
  Astar[i] <- step(x_pred_st[i]-cp[1])
  Cstar[i] <- 1+step(x_pred_st[i]-cp[2])
  mu_pred[i] <- alpha[Cstar[i]] + beta[Astar[i] + Cstar[i]]*(x_pred_st[i]-cp[Cstar[i]])
  }


}"
      myinitial <- function() {
        list(
          "alpha" = c(rnorm(2, 0, 3)),
          "beta" = c(rnorm(1, 0, 3), NA, rnorm(1, 0, 3)),
          "cp.temp" = c(runif(2, min(dat$x_st), max(dat$x_st)))
        )
      }
    }

    if (n_cp == 3) {
      run_model <- "

model{
####CP Regression model

  ###Data Loop
    for(j in 1:n_obs)
  {
  y[j] ~ dnorm(mu_y[j],tau[j])
  A[j] <- step(x[j]-cp[1])
  B[j] <- step(x[j]-cp[2])
  C[j] <- 1+step(x[j]-cp[3])
  mu_y[j] <- alpha[B[j] + C[j]] + beta[A[j] + B[j] + C[j]]*(x[j]-cp[B[j] + C[j]])
  tau[j] <- (y_err[j]^2 + sigma^2)^-1
  }

  ##Priors
  alpha[1] ~ dnorm(0.0,0.01)
  alpha[2] ~ dnorm(0.0,0.01)
  alpha[3] ~ dnorm(0.0,0.01)

  beta[1] ~ dnorm(0.0,0.01)
  beta[2] <- (alpha[2] - alpha[1])/(cp[2]-cp[1])
  beta[3] <- (alpha[3] - alpha[2])/(cp[3]-cp[2])
  beta[4]~dnorm(0.0,0.01)

  sigma ~ dt(0,4^-2,1)T(0,)

  cp.temp[1] ~ dunif(x_min,x_max)
  cp.temp[2] ~ dunif(x_min,x_max)
  cp.temp[3] ~ dunif(x_min,x_max)

  cp[1:3] <- sort(cp.temp)

  for(i in 1:n_pred)
  {
  Astar[i] <- step(x_pred_st[i]-cp[1])
  Bstar[i] <- step(x_pred_st[i]-cp[2])
  Cstar[i] <- 1+step(x_pred_st[i]-cp[3])
  mu_pred[i] <- alpha[Bstar[i] + Cstar[i]] + beta[Astar[i] + Bstar[i] + Cstar[i]]*(x_pred_st[i]-cp[Bstar[i] + Cstar[i]])
  }


}"
      myinitial <- function() {
        list(
          "alpha" = c(rnorm(3, 0, 3)),
          "beta" = c(rnorm(1, 0, 3), NA, NA, rnorm(1, 0, 3)),
          "cp.temp" = c(runif(3, min(dat$x_st), max(dat$x_st)))
        )
      }
    }
  }

  if (model == "cp" && EIV == TRUE) {
    if (n_cp == 1) {
      run_model <- "

model{
####CP Regression model

  ###Data Loop
    for(j in 1:n_obs)
  {
  y[j]~dnorm(mu_y[j],tau[j])
  x[j] ~ dnorm(mu_x[j],x_err[j]^-2)
  C[j] <- 1+step(mu_x[j]-cp)
  mu_x[j] ~ dnorm(0,0.5^-2)
  mu_y[j] <- alpha + beta[C[j]]*(mu_x[j]-cp)
  tau[j] <- (y_err[j]^2 + sigma^2)^-1

  }

  ##Priors
  alpha[1] ~ dnorm(0.0,0.01)

  beta[1]~dnorm(0.0,0.01)
  beta[2]~dnorm(0.0,0.01)

  sigma ~ dt(0,4^-2,1)T(0,)

  cp ~ dunif(x_min,x_max)

  for(i in 1:n_pred)
  {
  mu_pred[i] <-alpha + beta[Cstar[i]]*(x_pred_st[i]-cp)
  Cstar[i] <- 1+step(x_pred_st[i]-cp)
  }


}"
    }

    if (n_cp == 2) {
      run_model <- "

model{
####CP Regression model

  ###Data Loop
    for(j in 1:n_obs)
  {
  y[j] ~ dnorm(mu_y[j],tau[j])
  x[j] ~ dnorm(mu_x[j],x_err[j]^-2)
  A[j] <- step(mu_x[j]-cp[1])
  C[j] <- 1+step(mu_x[j]-cp[2])

  mu_x[j] ~ dnorm(0,0.5^-2)
  mu_y[j] <- alpha[C[j]] + beta[A[j] + C[j]]*(mu_x[j]-cp[C[j]])

  tau[j] <- (y_err[j]^2 + sigma^2)^-1

  }

  ##Priors
  alpha[1] ~ dnorm(0.0,0.01)
  alpha[2] ~ dnorm(0.0,0.01)

  beta[1]~dnorm(0.0,0.01)
  beta[2] <- (alpha[2] - alpha[1])/(cp[2]-cp[1])
  beta[3]~dnorm(0.0,0.01)

  sigma ~ dt(0,4^-2,1)T(0,)

  cp.temp[1] ~ dunif(x_min,x_max)
  cp.temp[2] ~ dunif(x_min,x_max)

  cp[1:2]<-sort(cp.temp)

  for(i in 1:n_pred)
  {
  Astar[i] <- step(x_pred_st[i]-cp[1])
  Cstar[i] <- 1+step(x_pred_st[i]-cp[2])
  mu_pred[i] <- alpha[Cstar[i]] + beta[Astar[i] + Cstar[i]]*(x_pred_st[i]-cp[Cstar[i]])
  }



}"

      myinitial <- function() {
        list(
          "alpha" = c(rnorm(2, 0, 3)),
          "beta" = c(rnorm(1, 0, 3), NA, rnorm(1, 0, 3)),
          "cp.temp" = c(runif(2, min(dat$x_st), max(dat$x_st)))
        )
      }
    }

    if (n_cp == 3) {
      run_model <- "

model{
####CP Regression model

  ###Data Loop
    for(j in 1:n_obs)
  {
  y[j] ~ dnorm(mu_y[j],tau[j])
  x[j] ~ dnorm(mu_x[j],x_err[j]^-2)

  A[j] <- step(mu_x[j]-cp[1])
  B[j] <- step(mu_x[j]-cp[2])
  C[j] <- 1+step(mu_x[j]-cp[3])

  mu_y[j] <- alpha[B[j] + C[j]] + beta[A[j] + B[j] + C[j]]*(mu_x[j]-cp[B[j] + C[j]])
  mu_x[j] ~ dnorm(0,0.5^-2)

  tau[j] <- (y_err[j]^2 + sigma^2)^-1
  }

  ##Priors
  alpha[1] ~ dnorm(0.0,0.01)
  alpha[2] ~ dnorm(0.0,0.01)
  alpha[3] ~ dnorm(0.0,0.01)

  beta[1]~dnorm(0.0,0.01)
  beta[2] <- (alpha[2] - alpha[1])/(cp[2]-cp[1])
  beta[3] <- (alpha[3] - alpha[2])/(cp[3]-cp[2])
  beta[4]~dnorm(0.0,0.01)

  sigma ~ dt(0,4^-2,1)T(0,)

  cp.temp[1] ~ dunif(x_min,x_max)
  cp.temp[2] ~ dunif(x_min,x_max)
  cp.temp[3] ~ dunif(x_min,x_max)

  cp[1:3]<-sort(cp.temp)

  for(i in 1:n_pred)
  {
  Astar[i] <- step(x_pred_st[i]-cp[1])
  Bstar[i] <- step(x_pred_st[i]-cp[2])
  Cstar[i] <- 1+step(x_pred_st[i]-cp[3])
  mu_pred[i] <- alpha[Bstar[i] + Cstar[i]] + beta[Astar[i] + Bstar[i] + Cstar[i]]*(x_pred_st[i]-cp[Bstar[i] + Cstar[i]])
  }


}"
      myinitial <- function() {
        list(
          "alpha" = c(rnorm(3, 0, 3)),
          "beta" = c(rnorm(1, 0, 3), NA, NA, rnorm(1, 0, 3)),
          "cp.temp" = c(runif(3, min(dat$x_st), max(dat$x_st)))
        )
      }
    }
  }

  if (model == "gp" && EIV == FALSE) {
    run_model <- "
model{

  gp ~ dmnorm(mu,Sigma.inv)
  Sigma.inv <- inverse(Sigma)
  mu_x <- x
  for(i in 1:n_obs)
  {
    mu[i] <- alpha
    Sigma[i,i] <- sigma_g^2 + 0.0001
    for(j in (i+1):n_obs) {
    Sigma[i,j] <- sigma_g^2*exp(-(phi^2)*((x[i]-x[j])^2))
    Sigma[j,i] <- Sigma[i,j]
    }

    y[i]~dnorm(gp[i],tau[i])
    tau[i] <- (y_err[i]^2 + sigma^2)^-1


  }


  sigma_g ~ dt(0,4^-2,1)T(0,)
  phi ~ dt(0,4^-2,4)T(0,)
  sigma ~ dt(0,4^-2,1)T(0,)
  alpha ~ dnorm(0,0.01)
}

"
  }

  if (model == "gp" && EIV == TRUE) {
    run_model <- "
model{

  gp ~ dmnorm(mu,Sigma.inv)
  Sigma.inv <- inverse(Sigma)

  for(i in 1:n_obs)
  {
    mu[i] <- alpha
    Sigma[i,i] <- sigma_g^2 + 0.0001
    for(j in (i+1):n_obs) {
    Sigma[i,j] <- sigma_g^2*exp(-(phi^2)*((mu_x[i]-mu_x[j])^2))
    Sigma[j,i] <- Sigma[i,j]
    }

    y[i]~dnorm(gp[i],tau[i])
    x[i] ~ dnorm(mu_x[i],x_err[i]^-2)
    mu_x[i] ~ dnorm(0,0.5^-2)
    tau[i] <- (y_err[i]^2 + sigma^2)^-1


  }


  sigma_g ~ dt(0,4^-2,1)T(0,)
  phi ~ dt(0,4^-2,4)T(0,)
  sigma ~ dt(0,4^-2,1)T(0,)
  alpha ~ dnorm(0,0.01)
}

"
  }

  if (model == "igp") {
    run_model <-
      "model
{
  for(i in 1:n_obs)
  {
    y[i]~dnorm(alpha + w.tilde.m[i],tau[i])
    x[i] ~ dnorm(mu_x[i],x_err[i]^-2)
    mu_x[i] ~ dnorm(0,0.5^-2)
    noisy_xerr[i] <- sqrt((beta^2)*(x_err[i]^2))
    tau[i] <- (y_err[i]^2 + sigma^2 + noisy_xerr[i]^2)^-1

  }


  ###Derivative process
  w.m~dmnorm(mu.w,K.inv)
  K.inv <- inverse((sigma_g^2)*K)
  K.w.inv<-inverse(K)

  for(i in 1:Ngrid)
  {
    mu.w[i] <- beta
    K[i,i]<-1+0.00001

    ######Exponential covariance for the derivative process
    for(j in (i+1):Ngrid)
    {
      K[i,j]<-(pow(phi,pow(Dist[i,j],1.99)))
      K[j,i]<-K[i,j]
    }
  }

  ###Expected value of the integrated Gaussian Process
  for(i in 1:Ngrid) {
    for(j in 1:n_obs) {
      K.gw[j,i]<-sum(pow(phi,quad1[j,i,])*quad2[j,i,])  #### Quadrature function
    } #End j loop
  } #End i loop

  w.tilde.m<-K.gw%*%K.w.inv%*%w.m

  # Priors
  sigma_g ~ dt(0,2^-2,1)T(0,)
  phi ~ dbeta(al,10)
  sigma ~ dt(0,2^-2,1)T(0,)
  alpha ~ dnorm(0,2^-2)
  beta ~ dnorm(0,2^-2)

}##End model
"
    igp_dat_list <- igp_data(dat)
  }


  ### The required data
  jags_data <- list(
    y = dat$y,
    y_err = dat$y_err,
    x = dat$x_st,
    x_err = dat$x_err_st,
    n_obs = nrow(dat),
    n_pred = 50,
    x_pred_st = seq(min(dat$x_st), max(dat$x_st), length.out = 50),
    x_pred = seq(min(dat$x), max(dat$x), length.out = 50),
    x_min = min(dat$x_st),
    x_max = max(dat$x_st),
    n_cp = n_cp,
    al = igp_smooth * 10 / (1 - igp_smooth)
  )


  ### Parameters to save
  if (model == "gp") {
    jags_pars <- c(
      "alpha",
      "phi",
      "sigma_g",
      "sigma",
      "mu_x"
    )
  }
  if (model == "slr") {
    jags_pars <- c(
      "mu_pred",
      "y_pred",
      "beta",
      "alpha",
      "sigma"
    )
  }
  if (model == "cp") {
    jags_pars <- c(
      "mu_pred",
      "y_pred",
      "beta",
      "alpha",
      "sigma",
      "cp"
    )
  }
  if (model == "igp") {
    jags_pars <- c(
      "phi",
      "sigma_g",
      "sigma",
      "w.m",
      "alpha",
      "beta"
    )
    jags_data <- c(igp_dat_list, jags_data)
  }



  ### Run the model
  mod <- suppressWarnings(R2jags::jags(
    data = jags_data,
    parameters.to.save = jags_pars,
    model.file = textConnection(run_model),
    n.iter = iter,
    n.burnin = burnin,
    n.thin = thin,
    inits = myinitial
  ))


  ### Create an object containing the posterior samples
  m <- mod$BUGSoutput$sims.matrix
  sample_draws <- tidybayes::tidy_draws(m)
  par_summary <- posterior::summarise_draws(sample_draws)
  if (sum(par_summary$rhat > 1.1, na.rm = TRUE) == 0) message("No convergence issues detected")
  if (sum(par_summary$rhat > 1.1, na.rm = TRUE) > 0) message("Convergence issues detected")


  return(list(
    sample_draws = sample_draws,
    model = model,
    EIV = EIV,
    dat = dat,
    jags_data = jags_data,
    scale_factor = scale_factor
  ))
}


#' Get parameter estimates and model estimates
#'
#' @param mod Fitted model object from \code{\link{run_mod}}
#'
#' @return List with model estimates (pred_summary) and parameter estimates (par_summary)
#' @export
#'
#' @examples
#' dat <- sim_slr(n_sim = 30)
#' mod <- run_mod(dat, model = "slr")
#' par_est(mod)
#'
par_est <- function(mod) {
  mu_pred <- .lower <- .upper <- x <- pred_y <- lwr_95 <- upr_95 <- alpha <- cp <- sigma_g <- phi <- sigma <- mu_x <- dat <- NULL

  sample_draws <- mod$sample_draws
  n_iter <- sample_draws$.iteration %>%
    unique() %>%
    length()
  if (n_iter > 1000) {
    sample_draws <- sample_draws %>% dplyr::slice_sample(n = 1000)
    n_iter <- 1000
  }
  jags_data <- mod$jags_data

  if (mod$model == "slr") {
    pred_summary <- sample_draws %>%
      tidyr::pivot_longer(`mu_pred[1]`:`mu_pred[50]`,
        names_to = "n",
        values_to = "mu_pred"
      ) %>%
      dplyr::mutate(n = rep(1:50, n_iter)) %>%
      dplyr::group_by(n) %>%
      tidybayes::median_qi(mu_pred) %>%
      dplyr::mutate(x = mod$jags_data$x_pred) %>%
      dplyr::mutate(
        pred_y = mu_pred,
        lwr_95 = .lower,
        upr_95 = .upper
      ) %>%
      dplyr::select(x, pred_y, lwr_95, upr_95)


    par_summary <- posterior::summarise_draws(sample_draws, ~quantile(.x,probs = c(0.025,0.975))) %>% dplyr::filter(variable %in% c("alpha", "beta", "sigma"))
  }

  if (mod$model == "cp") {

    ### Store results
    pred_summary <- sample_draws %>%
      tidyr::pivot_longer(`mu_pred[1]`:`mu_pred[50]`,
        names_to = "n",
        values_to = "mu_pred"
      ) %>%
      dplyr::mutate(n = rep(1:50, n_iter)) %>%
      dplyr::group_by(n) %>%
      tidybayes::median_qi(mu_pred) %>%
      dplyr::mutate(x = mod$jags_data$x_pred) %>%
      dplyr::mutate(
        pred_y = mu_pred,
        lwr_95 = .lower,
        upr_95 = .upper
      ) %>%
      dplyr::select(x, pred_y, lwr_95, upr_95)

    if (jags_data$n_cp == 1) {
      par_summary <- posterior::summarise_draws(sample_draws) %>%
        dplyr::filter(variable %in% c("alpha", "beta[1]", "beta[2]", "cp", "sigma"))
    }

    if (jags_data$n_cp == 2) {
      par_summary <- posterior::summarise_draws(sample_draws) %>%
        dplyr::filter(variable %in% c("alpha[1]", "alpha[2]", "beta[1]", "beta[2]", "beta[3]", "cp[1]", "cp[2]", "sigma"))
    }

    if (jags_data$n_cp == 3) {
      par_summary <- posterior::summarise_draws(sample_draws) %>%
        dplyr::filter(variable %in% c("alpha[1]", "alpha[2]", "alpha[3]", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "cp[1]", "cp[2]", "cp[3]", "sigma"))
    }
  }

  if (mod$model == "gp") {
    n_pred <- 50
    n_obs <- jags_data$n_obs
    index <- 1:n_obs
    x_star <- seq(min(jags_data$x), max(jags_data$x), length.out = n_pred)
    par_est <- posterior::summarise_draws(sample_draws)

    # posterior estimate for pars
    sigma <- par_est %>%
      dplyr::filter(variable == "sigma") %>%
      dplyr::pull(median)
    alpha <- par_est %>%
      dplyr::filter(variable == "alpha") %>%
      dplyr::pull(median)
    phi <- par_est %>%
      dplyr::filter(variable == "phi") %>%
      dplyr::pull(median)
    sigma_g <- par_est %>%
      dplyr::filter(variable == "sigma_g") %>%
      dplyr::pull(median)
    mu_x <- par_est %>%
      dplyr::filter(!variable %in% c("sigma", "alpha", "phi", "sigma_g", "deviance")) %>%
      dplyr::pull(median)

    ### Predicitive distribution for the GP
    Sigma <- sigma * 2 * diag(n_obs) + sigma_g^2 * exp(-(phi^2) * fields::rdist(mu_x, mu_x)^2)
    Sigma_star <- sigma_g^2 * exp(-(phi^2) * fields::rdist(x_star, mu_x)^2)
    Sigma_star_star <- sigma_g^2 * exp(-(phi^2) * fields::rdist(x_star, x_star)^2)


    pred_mean <- Sigma_star %*% solve(Sigma, mod$dat$y)
    pred_var <- Sigma_star_star - Sigma_star %*% solve(Sigma, t(Sigma_star))

    temp <- mvtnorm::rmvnorm(n_iter, pred_mean, pred_var)
    deriv <- matrix(NA, nrow = n_iter / 2, ncol = length(jags_data$x_pred_st))
    newD <- jags_data$x_pred_st

    for (i in 1:(n_iter / 2))
    {
      dat <- tibble::tibble(x = mod$jags_data$x_pred_st, y = temp[i, ])
      gam_mod <- mgcv::gamm(y ~ s(x, k = 30), data = dat)
      gam_mod <- gam_mod$gam
      X0 <- mgcv::predict.gam(gam_mod, data.frame(x = newD), type = "lpmatrix")
      newD <- newD + 1e-7
      X1 <- mgcv::predict.gam(gam_mod, data.frame(x = newD), type = "lpmatrix")
      Xp <- (X1 - X0) / 1e-7
      Xp.r <- NROW(Xp)
      Xp.c <- NCOL(Xp)
      # number of smooth terms
      Xi <- Xp * 0
      want <- grep("x", colnames(X1))
      Xi[, want] <- Xp[, want]
      df <- Xi %*% coef(gam_mod)
      deriv[i, ] <- c(df)
    }

    ### Store results
    pred_summary <- tibble::tibble(
      x = x_star * mod$scale_factor,
      pred_y = c(pred_mean),
      lwr_95 = pred_y - 1.96 * sqrt(diag(pred_var)),
      upr_95 = pred_y + 1.96 * sqrt(diag(pred_var)),
      rate_y = apply(deriv, 2, median),
      rate_lwr_95 = apply(deriv, 2, quantile, probs = 0.025),
      rate_upr_95 = apply(deriv, 2, quantile, probs = 0.975)
    )

    par_summary <- posterior::summarise_draws(sample_draws) %>%
      dplyr::filter(variable %in% c("alpha", "phi", "sigma_g", "sigma"))
  }

  if (mod$model == "igp") {

    # Get predictions on a grid of x values.
    Ngrid <- jags_data$Ngrid
    xgrid <- jags_data$xstar
    xstar <- jags_data$xstar

    Dist <- jags_data$Dist

    # Set up the matrix that will contain the estimates
    pred <- matrix(NA, ncol = Ngrid, nrow = n_iter)
    K.gw <- K <- K.w.inv <- array(NA, c(n_iter, Ngrid, Ngrid))

    ######## Initialize quadrature for the integration########
    L <- 30 ## this sets the precision of the integration quadrature (higher is better but more computationally expensive)
    index <- 1:L
    cosfunc <- cos(((2 * index - 1) * pi) / (2 * L))

    quad1 <- array(dim = c(nrow = Ngrid, ncol = Ngrid, L))
    quad2 <- array(dim = c(nrow = Ngrid, ncol = Ngrid, L))

    for (j in 1:Ngrid)
    {
      for (k in 1:Ngrid)
      {
        quad1[k, j, ] <- abs((xgrid[k] * cosfunc / 2) + (xgrid[k] / 2) - xstar[j])^1.99
        quad2[k, j, ] <- ((xgrid[k] / 2) * (pi / L)) * (sqrt(1 - cosfunc^2))
      }
    }

    # Get posterior samples of rates
    w.ms <- as.matrix(sample_draws %>% dplyr::select(`w.m[1]`:`w.m[50]`))

    # Get estimates
    for (i in 1:n_iter) {
      for (k in 1:Ngrid) {
        for (j in 1:Ngrid) {
          K.gw[i, j, k] <- sum((sample_draws$phi[i]^quad1[j, k, ]) * quad2[j, k, ]) #### Quadrature function
        } # End j loop
      } # End k loop

      K[i, , ] <- sample_draws$phi[i]^(Dist^1.99)
      K.w.inv[i, , ] <- solve(K[i, , ])
      pred[i, ] <- sample_draws$alpha[i] + K.gw[i, , ] %*% K.w.inv[i, , ] %*% w.ms[i, ]
    } # End i loop

    ### Store results
    pred_summary <- tibble::tibble(
      x = seq(min(mod$dat$x), max(mod$dat$x), length.out = 50),
      pred_y = apply(pred, 2, stats::median),
      lwr_95 = apply(pred, 2, stats::quantile, probs = 0.025),
      upr_95 = apply(pred, 2, stats::quantile, probs = 0.975),
      rate_y = apply(w.ms, 2, stats::median),
      rate_lwr_95 = apply(w.ms, 2, stats::quantile, probs = 0.025),
      rate_upr_95 = apply(w.ms, 2, stats::quantile, probs = 0.975)
    )

    par_summary <- posterior::summarise_draws(sample_draws) %>%
      dplyr::filter(variable %in% c("phi", "sigma_g", "sigma"))
  }

  return(list(
    pred_summary = pred_summary,
    par_summary = par_summary
  ))
}
