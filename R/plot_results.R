#' Plot of data plus model estimates (and credible intervals)
#' @param mod1 The fitted model
#' @param mod2 An optional second fitted model for comparison
#' @param add_truth Logical argument to add the "True" data generating process to the plot. This should only be set to True when using simulated data from sim_dat
#' @return Plot
#' @export
#' @examples
#' dat <- sim_slr(n_sim = 30)
#' mod <- run_mod(dat,model = "model_eiv_slr")
#' plot_res(mod)
plot_res <- function(mod1,
                     mod2 = NULL,
                     add_truth = FALSE) {

  x <- pred_y <- lwr_95 <- upr_95 <- rate_y <- rate_lwr_95 <- rate_upr_95 <- y <- obs_index <- true_x <- true_y <- model_label <- NULL
  ### format data for plotting with errors
  data_to_plot <- dat_to_plot(mod1$dat)

  ## get model based estimates
  res <- par_est(mod1)
  pred_res <- res$pred_res
  pred_res$model_label <- mod1$model

  ### same set up if second model is included
  if (!is.null(mod2)) {
    pred_res1 <- pred_res
    res2 <- par_est(mod2)
    pred_res2 <- res2$pred_res
    pred_res1$model_label <- mod1$model
    pred_res2$model_label <- mod2$model
    pred_res <- dplyr::full_join(pred_res1, pred_res2)
  }

  ### plot the results
  p <- ggplot2::ggplot(pred_res, ggplot2::aes(x = x, y = pred_y)) +
    ggplot2::geom_line(ggplot2::aes(colour = model_label)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr_95, ymax = upr_95, fill = model_label), alpha = 0.4) +
    ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, group = obs_index), data = data_to_plot, alpha = 0.2) +
    ggplot2::ylab("Y") +
    ggplot2::xlab("X") +
    ggplot2::labs(colour = "", fill = "95% UI") +
    ggplot2::theme_classic()

  if(mod1$model == "model_eiv_igp")
  {
  p_rate <- ggplot2::ggplot(pred_res, ggplot2::aes(x = x, y = rate_y)) +
      ggplot2::geom_line(ggplot2::aes(colour = model_label)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = rate_lwr_95, ymax = rate_upr_95, fill = model_label), alpha = 0.4) +
      ggplot2::ylab("Rate") +
      ggplot2::xlab("X") +
      ggplot2::labs(colour = "", fill = "95% UI") +
      ggplot2::theme_classic()

  p <- ggpubr::ggarrange(p,p_rate,ncol = 1, common.legend = TRUE,align = "h")
  }
  ### add true line if indicated
  if (add_truth) {
    return(p + ggplot2::geom_line(data = mod1$dat, ggplot2::aes(x = true_x, y = true_y, colour = "Truth"), size = 1.5) +
             ggplot2::labs(colour = ""))
  }

  if (!add_truth) {
    return(p)
  }
}
