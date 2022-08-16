#' Plot data with measurement uncertainty (in y and x variables)
#'
#' @param dat Input data with columns x,x_err,y,y_err
#' @param add_truth Logical argument to add the "True" data generating process to the plot. This should only be set to True when using simulated data from sim_slr, sim_cp or sim_gp.
#'
#' @return Plot of data with measurement errors
#' @export
#'
#' @examples
#' dat <- sim_slr(n_sim = 30)
#' plot_dat(dat)
plot_dat <- function(dat,
                     BP_scale = FALSE,
                     add_truth = FALSE) {

  x <- y <- obs_index <- true_x <- true_y <- NULL
  data_to_plot <- dat_to_plot(dat)
  p <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_polygon(ggplot2::aes(group = obs_index, fill = "1-sigma error"), alpha = 0.3) +
    ggplot2::geom_point(data = dat, ggplot2::aes(x = x, y = y), alpha = 0.6, pch = 1) +
    ggplot2::xlab("X") +
    ggplot2::ylab("Y") +
    ggplot2::scale_fill_manual(values = "gray", name = "") +
    ggplot2::theme_classic()


  if(BP_scale) p <- p + scale_x_reverse()

  if (add_truth) {
    return(p + ggplot2::geom_line(data = dat, ggplot2::aes(x = true_x, y = true_y, colour = "Truth"), size = 1.5) +
             ggplot2::labs(colour = ""))
  }

  if (!add_truth) {
    return(p)
  }
}


dat_to_plot <- function(raw_dat) {

  y <- y_err <- x <- x_err <- y_1_lwr <- x_4_upr <- bounds <- NULL

  dat <- raw_dat %>%
    dplyr::mutate(
      y_1_lwr = y - y_err,
      y_2_lwr = y - y_err,
      y_3_upr = y + y_err,
      y_4_upr = y + y_err,
      x_1_upr = x + x_err,
      x_2_lwr = x - x_err,
      x_3_lwr = x - x_err,
      x_4_upr = x + x_err
    )


  get_bounds <- dat %>%
    dplyr::select(y_1_lwr:x_4_upr) %>%
    dplyr::mutate(obs_index = 1:dplyr::n()) %>%
    tidyr::pivot_longer(
      cols = y_1_lwr:x_4_upr,
      names_to = "bounds",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      bounds = replace(bounds, bounds %in% c("y_1_lwr", "y_2_lwr", "y_3_upr", "y_4_upr"), "y"),
      bounds = replace(bounds, bounds %in% c("x_1_upr", "x_2_lwr", "x_3_lwr", "x_4_upr"), "x")
    )

  x_bounds <- get_bounds %>%
    dplyr::filter(bounds == "x")

  y_bounds <- get_bounds %>%
    dplyr::filter(bounds == "y")

  data_to_plot <- tibble::tibble(
    obs_index = x_bounds$obs_index,
    x = x_bounds$value,
    y = y_bounds$value
  )

  return(data_to_plot)
}
