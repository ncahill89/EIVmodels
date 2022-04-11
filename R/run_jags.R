#' Fit errors-in-variables models with JAGS
#'
#' @param dat Input data with columns x,x_err,y,y_err
#' @param model The model to run. Choose from model_reg, model_eiv_reg, model_cp, model_eiv_cp, model_gp, model_eiv_gp, model_eiv_igp. Defaults to model_eiv_reg.
#' @param iter MCMC iterations
#' @param burnin MCMC burnin
#' @param thin MCMC thinning
#' @param scale_factor value to divide the predictor (x) by to change the scale. 1 will result in n0 change. 1000 is recommended if x is years.
#'
#' @return List with data, JAGS input data, model output and the name of the model file.
#' @export
#'
#' @examples
#' dat <- sim_reg(n_sim = 30)
#' mod <- run_mod(dat, model = "model_eiv_reg")
run_mod <- function(dat,
                    model = "model_eiv_reg",
                    iter = 15000,
                    burnin = 5000,
                    thin = 5,
                    scale_factor = 1) {

  # Simple Linear Regression Model ------------------------------------------
  x<-x_err<- NULL

  dat <- dat %>% dplyr::mutate(x_st = x/scale_factor,
                               x_err_st = x_err/scale_factor)

  if (model == "model_reg") {
    run_model <-
      "model{
  ## Data model loop
  for(j in 1:n_obs)
  {
  y[j]~dnorm(mu_y[j],(sigma + y_err[j])^-2)
  mu_y[j] <- alpha + beta*(x[j])
  } # end j loop

 ## Priors
 alpha ~ dnorm(0.0,0.01)
 beta ~ dnorm(0.0,0.01)
 sigma ~ dt(0,4^-2,1)T(0,)

 for(i in 1:n_pred) {mu_pred[i] = alpha + beta*x_pred_st[i]}

}# end
"
  }

  if (model == "model_eiv_reg") {
    run_model <-
      "model{
  ## Data model loop
  for(j in 1:n_obs)
  {
  y[j]~dnorm(mu_y[j],(sigma + y_err[j])^-2)
  x[j] ~ dnorm(mu_x[j],x_err[j]^-2)
  mu_x[j] ~ dnorm(0,0.5^-2)
  mu_y[j] <- alpha + beta*(mu_x[j])
  } # end j loop

 ## Priors
 alpha ~ dnorm(0.0,0.01)
 beta ~ dnorm(0.0,0.01)
 sigma ~ dt(0,4^-2,1)T(0,)

 for(i in 1:n_pred) {mu_pred[i] = alpha + beta*x_pred_st[i]}

}# end
"
  }

  if (model == "model_cp1") {
    run_model <- "

model{
####CP Regression model

  ###Data Loop
    for(j in 1:n_obs)
  {
  y[j]~dnorm(mu_y[j],(sigma + y_err[j])^-2)
  C[j] <- 1+step(x[j]-cp)
  mu_y[j] <- alpha + beta[C[j]]*(x[j]-cp)
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

  if (model == "model_eiv_cp1") {
    run_model <- "

model{
####CP Regression model

  ###Data Loop
    for(j in 1:n_obs)
  {
  y[j]~dnorm(mu_y[j],(sigma + y_err[j])^-2)
  x[j] ~ dnorm(mu_x[j],x_err[j]^-2)
  C[j] <- 1+step(mu_x[j]-cp)
  mu_x[j] ~ dnorm(0,0.5^-2)
  mu_y[j] <- alpha + beta[C[j]]*(mu_x[j]-cp)
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

  if (model == "model_gp") {
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

    y[i]~dnorm(gp[i],(sigma + y_err[i])^-2)

  }


  sigma_g ~ dt(0,4^-2,1)T(0,)
  phi ~ dt(0,4^-2,4)T(0,)
  sigma ~ dt(0,4^-2,1)T(0,)
  alpha ~ dnorm(0,0.01)
}

"
  }

  if (model == "model_eiv_gp") {
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

    y[i]~dnorm(gp[i],(sigma + y_err[i])^-2)
    x[i] ~ dnorm(mu_x[i],x_err[i]^-2)
    mu_x[i] ~ dnorm(0,0.5^-2)

  }


  sigma_g ~ dt(0,4^-2,1)T(0,)
  phi ~ dt(0,4^-2,4)T(0,)
  sigma ~ dt(0,4^-2,1)T(0,)
  alpha ~ dnorm(0,0.01)
}

"
  }

  if (model == "model_eiv_igp") {
  run_model <-
"model
{
  for(i in 1:n_obs)
  {
    y[i]~dnorm(alpha + w.tilde.m[i],(sigma + y_err[i] + noisy_xerr[i])^-2)
    x[i] ~ dnorm(mu_x[i],x_err[i]^-2)
    mu_x[i] ~ dnorm(0,0.5^-2)
    noisy_xerr[i] <- sqrt((beta^2)*(x_err[i]^2))
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
  phi ~ dunif(0,1)
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
    x_max = max(dat$x_st)
    )


  ### Parameters to save
  if (model == "model_gp" | model == "model_eiv_gp") {
    jags_pars <- c(
      "alpha",
      "phi",
      "sigma_g",
      "sigma",
      "mu_x"
    )
  }
  if (model == "model_reg" | model == "model_eiv_reg") {
    jags_pars <- c(
      "mu_pred",
      "y_pred",
      "beta",
      "alpha",
      "sigma"
    )
  }
  if (model == "model_cp1" | model == "model_eiv_cp1") {
    jags_pars <- c(
      "mu_pred",
      "y_pred",
      "beta",
      "alpha",
      "sigma",
      "cp"
    )
  }
  if (model == "model_eiv_igp") {

    jags_pars <- c(
      "phi",
      "sigma_g",
      "sigma",
      "w.m",
      "noisy_xerr",
      "alpha",
      "beta"
    )
  jags_data <- c(igp_dat_list,jags_data)
  }



  ### Run the model
  mod <- suppressWarnings(R2jags::jags(
    data = jags_data,
    parameters.to.save = jags_pars,
    model.file = textConnection(run_model),
    n.iter = iter,
    n.burnin = burnin,
    n.thin = thin
  ))


  ### Create an object containing the posterior samples
  m <- mod$BUGSoutput$sims.matrix


  return(list(
    m = m,
    model = model,
    dat = dat,
    sims_list = mod$BUGSoutput$sims.list,
    jags_data = jags_data,
    scale_factor = scale_factor
  ))
  }


#' Get parameter estimates and model estimates
#'
#' @param mod Fitted model object
#'
#' @return List with model estimates (pred_res) and parameter estimates (par_summary)
#' @export
#'
#' @examples
#' dat <- sim_reg(n_sim = 30)
#' mod <- run_mod(dat, model = "model_eiv_reg")
#' par_est(mod)

par_est <- function(mod) {

  mu_pred <- .lower <- .upper <- x <- pred_y <- lwr_95 <- upr_95 <- alpha <- cp <- sigma_g <- phi <- sigma <- mu_x <- dat <- NULL

  if (mod$model == "model_reg" | mod$model == "model_eiv_reg") {
    pred_res <- mod$m %>%
      tidybayes::spread_draws(mu_pred[1:mod$jags_data$n_pred]) %>%
      tidybayes::median_qi(mu_pred) %>%
      dplyr::mutate(x = mod$jags_data$x_pred) %>%
      dplyr::mutate(
        pred_y = mu_pred,
        lwr_95 = .lower,
        upr_95 = .upper
      ) %>%
      dplyr::select(x, pred_y, lwr_95, upr_95)


    par_dat <- mod$m %>%
      tidybayes::gather_draws(alpha, beta,sigma) %>%
      tidybayes::median_qi()
  }

  if (mod$model == "model_cp1" | mod$model == "model_eiv_cp1") {

    ### Store results
    pred_res <- mod$m %>%
      tidybayes::spread_draws(mu_pred[1:mod$jags_data$n_pred]) %>%
      tidybayes::median_qi(mu_pred) %>%
      dplyr::mutate(x = mod$jags_data$x_pred) %>%
      dplyr::mutate(
        pred_y = mu_pred,
        lwr_95 = .lower,
        upr_95 = .upper
      ) %>%
      dplyr::select(x, pred_y, lwr_95, upr_95)

    par_dat <- mod$m %>%
      tidybayes::gather_draws(alpha, beta[1:2], cp,sigma) %>%
      tidybayes::median_qi()
  }

  if (mod$model == "model_gp" | mod$model == "model_eiv_gp") {
    m <- mod$m
    sims_list <- mod$sims_list
    jags_data <- mod$jags_data

    n_pred <- 50
    n_obs <- jags_data$n_obs
    index <- 1:n_obs
    x_star <- seq(min(jags_data$x), max(jags_data$x), length.out = n_pred)
    par_dat <- m %>% tidybayes::spread_draws(sigma_g, phi, sigma, mu_x[index]) # posterior samples
    par_est <- par_dat %>%
      tidybayes::mean_qi(sigma_g, phi, sigma, mu_x) # posterior estimate for pars

    ### Predicitive distribution for the GP

    Sigma <- unique(par_est$sigma) * 2 * diag(n_obs) + unique(par_est$sigma_g)^2 * exp(-(unique(par_est$phi)^2) * fields::rdist(par_est$mu_x, par_est$mu_x)^2)
    Sigma_star <- unique(par_est$sigma_g)^2 * exp(-(unique(par_est$phi)^2) * fields::rdist(x_star, par_est$mu_x)^2)
    Sigma_star_star <- unique(par_est$sigma_g)^2 * exp(-(unique(par_est$phi)^2) * fields::rdist(x_star, x_star)^2)


    pred_mean <- Sigma_star %*% solve(Sigma, mod$dat$y)
    pred_var <- Sigma_star_star - Sigma_star %*% solve(Sigma, t(Sigma_star))

    ### Store results
    pred_res <- tibble::tibble(
      x = x_star*mod$scale_factor,
      pred_y = c(pred_mean),
      lwr_95 = pred_y - 1.96 * sqrt(diag(pred_var)),
      upr_95 = pred_y + 1.96 * sqrt(diag(pred_var))
    )

    par_dat <- mod$m %>%
      tidybayes::gather_draws(sigma,sigma_g, phi) %>%
      tidybayes::median_qi()
  }

  if (mod$model == "model_eiv_igp") {

    m <- mod$m
    sims_list <- mod$sims_list
    jags_data <- mod$jags_data
    n_iter <- min(dim(sims_list$w.m)[1],500)


    #Get predictions on a grid of x values.
    Ngrid <- jags_data$Ngrid
    xgrid <- jags_data$xstar
    xstar <- jags_data$xstar

    Dist <- jags_data$Dist

    #Set up the matrix that will contain the estimates
    pred <- matrix(NA,ncol=Ngrid,nrow=n_iter)
    K.gw<-K<-K.w.inv<-array(NA,c(n_iter, Ngrid, Ngrid))

    ########Initialize quadrature for the integration########
    L=30    ## this sets the precision of the integration quadrature (higher is better but more computationally expensive)
    index=1:L
    cosfunc=cos(((2*index-1)*pi)/(2*L))

    quad1=array(dim=c(nrow=Ngrid,ncol=Ngrid,L))
    quad2=array(dim=c(nrow=Ngrid,ncol=Ngrid,L))

    for(j in 1:Ngrid)
    {
      for(k in 1:Ngrid)
      {
        quad1[k,j,]=abs((xgrid[k]*cosfunc/2)+(xgrid[k]/2)-xstar[j])^1.99
        quad2[k,j,]=((xgrid[k]/2)*(pi/L))*(sqrt(1-cosfunc^2))
      }
    }

    #Get posterior samples of rates

    w.ms <- sims_list$w.m[1:n_iter,]


    #Get estimates
    for(i in 1:n_iter) {
      for(k in 1:Ngrid) {
        for(j in 1:Ngrid) {
          K.gw[i,j,k]<-sum((sims_list$phi[i]^quad1[j,k,])*quad2[j,k,])  #### Quadrature function
        } #End j loop
      } #End k loop

      K[i,,]<-sims_list$phi[i]^(Dist^1.99)
      K.w.inv[i,,]<-solve(K[i,,])
      pred[i,] <- sims_list$alpha[i] + K.gw[i,,]%*%K.w.inv[i,,]%*%w.ms[i,]
    } #End i loop


    ### Store results
    pred_res <- tibble::tibble(
      x = seq(min(mod$dat$x), max(mod$dat$x), length.out = 50),
      pred_y = apply(pred, 2, stats::median),
      lwr_95 = apply(pred, 2, stats::quantile, probs = 0.025),
      upr_95 = apply(pred, 2, stats::quantile, probs = 0.975))

    par_dat <- mod$m %>%
      tidybayes::gather_draws(sigma, sigma_g, phi) %>%
      tidybayes::median_qi()

  }

  return(list(
    pred_res = pred_res,
    par_summary = par_dat
  ))
}