#' Plot of data plus model estimates (and credible intervals)
#' @param mod1 The fitted model from \code{\link{run_mod}}
#' @param mod2 An optional second fitted model from \code{\link{run_mod}} for comparison
#' @param add_truth Logical argument to add the "True" data generating process to the plot. This should only be set to True when using simulated data from sim_dat
#' @return p
#' @export
#' @examples
#' dat <- sim_slr(n_sim = 30)
#' mod <- run_mod(dat, model = "slr")
#' plot_res(mod)
plot_res <- function(mod1,
                     mod2 = NULL,
                     add_truth = FALSE) {

  x <- pred_y <- lwr_95 <- upr_95 <- rate_y <- rate_lwr_95 <- rate_upr_95 <- y <- obs_index <- true_x <- true_y <- model_label <- NULL
  ### format data for plotting with errors
  data_to_plot <- dat_to_plot(mod1$dat)

  ## get model based estimates
  res <- par_est(mod1)
  pred_res <- res$pred_summary
  pred_res$model_label <- ifelse(mod1$EIV,paste0("eiv-",mod1$model),mod1$model)
  BP_scale <- mod1$BP_scale

  ### same set up if second model is included
  if (!is.null(mod2)) {
    pred_res1 <- pred_res
    res2 <- par_est(mod2)
    pred_res2 <- res2$pred_summary
    BP_scale <- mod2$BP_scale
    pred_res1$model_label <- ifelse(mod1$EIV, paste0("eiv-",mod1$model),mod1$model)
    pred_res2$model_label <- ifelse(mod2$EIV, paste0("eiv-",mod2$model),mod2$model)
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

  if(BP_scale) p <- p + ggplot2::scale_x_reverse()

  if(ncol(pred_res) == 8)
  {
  p_rate <- ggplot2::ggplot(pred_res %>% tidyr::drop_na(), ggplot2::aes(x = x, y = rate_y)) +
      ggplot2::geom_line(ggplot2::aes(colour = model_label)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = rate_lwr_95, ymax = rate_upr_95, fill = model_label), alpha = 0.4) +
      ggplot2::ylab("Rate") +
      ggplot2::xlab("X") +
      ggplot2::labs(colour = "", fill = "95% UI") +
      ggplot2::theme_classic()

  if(BP_scale) p_rate <- p_rate + ggplot2::scale_x_reverse()

  }
  ### add true line if indicated
  if(ncol(pred_res) == 8)
  {
  if (add_truth) {
    return(list(p = p + ggplot2::geom_line(data = mod1$dat, ggplot2::aes(x = true_x, y = true_y, colour = "Truth"), size = 1.5) +
             ggplot2::labs(colour = ""),
             p_rate = p_rate))
  }

  if (!add_truth) {
    return(list(p = p,
                p_rate = p_rate))
  }
  }

  if(mod1$model != "igp")
  {
    if (add_truth) {
      return(p + ggplot2::geom_line(data = mod1$dat, ggplot2::aes(x = true_x, y = true_y, colour = "Truth"), size = 1.5) +
                    ggplot2::labs(colour = ""))
    }

    if (!add_truth) {
      return(p)
    }
  }
}
