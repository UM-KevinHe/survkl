#' Plot Model Performance vs Eta for `coxkl`
#'
#' @description
#' Plots model performance across the \code{eta} sequence. Performance is either
#' loss (\code{-2} times partial log-likelihood) or concordance index (C-index).
#' If no test data are provided, the curve is computed on the training data stored
#' in \code{x$data}.
#'
#' @param x A fitted model object of class \code{"coxkl"}.
#' @param test_z Optional numeric matrix of test covariates.
#' @param test_time Optional numeric vector of test survival times.
#' @param test_delta Optional numeric vector of test event indicators.
#' @param test_stratum Optional vector of test stratum membership.
#' @param criteria Character string: \code{"loss"} or \code{"CIndex"}.
#' @param ... Additional arguments (ignored).
#' 
#' @details
#' When \code{criteria = "loss"} and no test data are supplied, the plotted values are
#' \code{(-2 * x$likelihood) / n}, where \code{n} is the number of rows in the
#' (training) data. When test data are provided, performance is computed via
#' \code{test_eval(..., criteria = "loss")} and divided by the test sample size.
#' For \code{criteria = "CIndex"}, performance is computed via
#' \code{test_eval(..., criteria = "CIndex")} on the chosen dataset. The plot adds a
#' dotted horizontal reference line at the value corresponding to \code{eta = 0}
#' (closest point on the \code{eta} grid).
#'
#' @return A \code{ggplot} object showing the performance curve.
#' 
#' @examples
#' data(ExampleData_lowdim)
#' 
#' train_dat_lowdim  <- ExampleData_lowdim$train
#' test_dat_lowdim   <- ExampleData_lowdim$test
#' beta_external_good_lowdim <- ExampleData_lowdim$beta_external_good
#' eta_grid <- generate_eta(method = "exponential", n = 100, max_eta = 30)
#' 
#' model <- coxkl(z = train_dat_lowdim$z,
#'                delta = train_dat_lowdim$status,
#'                time = train_dat_lowdim$time,
#'                stratum = train_dat_lowdim$stratum,
#'                beta = beta_external_good_lowdim,
#'                etas = eta_grid)
#' plot(model,
#'      test_z = test_dat_lowdim$z, 
#'      test_time = test_dat_lowdim$time, 
#'      test_delta = test_dat_lowdim$status, 
#'      test_stratum = test_dat_lowdim$stratum, 
#'      criteria = "loss")
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_segment labs coord_cartesian theme_minimal theme element_blank element_line element_text scale_color_manual scale_linetype_manual guides guide_legend
#' @importFrom grid unit
#' @importFrom rlang .data
#' @exportS3Method plot coxkl
plot.coxkl <- function(x, test_z = NULL, test_time = NULL, test_delta = NULL,
                       test_stratum = NULL, criteria = c("loss", "CIndex"), ...) {
  criteria <- match.arg(criteria)
  if (!inherits(x, "coxkl")) stop("'x' must be of class 'coxkl'.", call. = FALSE)
  
  etas          <- x$eta
  beta_mat      <- x$beta
  z_train       <- x$data$z
  time_train    <- x$data$time
  delta_train   <- x$data$delta
  stratum_train <- x$data$stratum
  
  using_train <- is.null(test_z) && is.null(test_time) && is.null(test_delta) && is.null(test_stratum)
  if (using_train) {
    test_z       <- z_train
    test_time    <- time_train
    test_delta   <- delta_train
    test_stratum <- stratum_train
  }
  
  n_eval <- nrow(as.matrix(test_z))
  
  if (criteria == "loss") {
    if (using_train) {
      metrics <- (-2 * x$likelihood) / n_eval
    } else {
      raw_metrics <- sapply(seq_along(etas), function(i)
        test_eval(
          test_z       = test_z,
          test_delta   = test_delta,
          test_time    = test_time,
          test_stratum = test_stratum,
          betahat      = beta_mat[, i],
          criteria     = "loss"
        )
      )
      metrics <- as.numeric(raw_metrics) / n_eval
    }
  } else {
    metrics <- sapply(seq_along(etas), function(i)
      test_eval(
        test_z       = test_z,
        test_delta   = test_delta,
        test_time    = test_time,
        test_stratum = test_stratum,
        betahat      = beta_mat[, i],
        criteria     = "CIndex"
      )
    )
  }
  
  idx0 <- which.min(abs(etas - 0))
  x0   <- etas[idx0]
  y0   <- as.numeric(metrics[idx0])
  xmax <- max(etas, na.rm = TRUE)
  
  df   <- data.frame(eta = etas, metric = as.numeric(metrics))
  ylab <- if (criteria == "CIndex") "C Index" else "Loss"
  
  p_main <- ggplot2::ggplot(df, ggplot2::aes(x = .data$eta, y = .data$metric)) +
    ggplot2::geom_line(size = 1, color = "#7570B3") +
    ggplot2::geom_point(size = 2, color = "#7570B3") +
    ggplot2::geom_segment(
      data = data.frame(x0 = x0, y0 = y0, xmax = xmax),
      ggplot2::aes(x = .data$x0, xend = .data$xmax, y = .data$y0, yend = .data$y0),
      inherit.aes = FALSE, color = "#1B9E77", linetype = "dotted", linewidth = 1
    ) +
    ggplot2::geom_point(
      data = data.frame(x0 = x0, y0 = y0),
      ggplot2::aes(x = .data$x0, y = .data$y0),
      inherit.aes = FALSE, color = "#1B9E77", shape = 16, size = 3
    ) +
    ggplot2::labs(x = expression(eta), y = ylab) +
    theme_biometrics() +
    ggplot2::coord_cartesian(
      ylim = c(min(c(df$metric, y0), na.rm = TRUE) * 0.995,
               max(c(df$metric, y0), na.rm = TRUE) * 1.005)
    )
  
  legend_df <- data.frame(
    x = rep(c(0, 1), 2),
    y = rep(1, 4),
    Method = factor(rep(c("survkl", "Internal"), each = 2),
                    levels = c("survkl", "Internal"))
  )
  
  internal_df <- legend_df[legend_df$Method == "Internal", , drop = FALSE]
  
  g_legend <- ggplot2::ggplot(legend_df, ggplot2::aes(x = .data$x, y = .data$y, 
                                                      color = .data$Method, linetype = .data$Method)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(
      data = internal_df,
      ggplot2::aes(x = 0.5, y = 1, color = .data$Method),
      inherit.aes = FALSE, shape = 16, size = 2.4
    ) +
    ggplot2::scale_color_manual(values = c("survkl" = "#7570B3", "Internal" = "#1B9E77")) +
    ggplot2::scale_linetype_manual(values = c("survkl" = "solid", "Internal" = "dotted")) +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(legend.position = "top",
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 14),
                   legend.key.width = grid::unit(1.2, "lines"),
                   legend.key.height = grid::unit(0.6, "lines")) +
    ggplot2::guides(color = ggplot2::guide_legend(keywidth = 1.2, keyheight = 0.4, title = NULL),
                    linetype = ggplot2::guide_legend(keywidth = 1.2, keyheight = 0.4, title = NULL))
  
  cowplot::plot_grid(cowplot::get_legend(g_legend), p_main, ncol = 1, rel_heights = c(0.08, 1))
}



#' Plot Model Performance vs Lambda for `coxkl_ridge`
#'
#' @description
#' Plots model performance across the \code{lambda} sequence. Performance is
#' loss (\code{-2} times partial log-likelihood) or concordance index (C-index).
#' If no test data are provided, the curve uses the training data stored in \code{x$data}.
#'
#' @param x A fitted model object of class \code{"coxkl_ridge"}.
#' @param test_z Optional numeric matrix of test covariates.
#' @param test_time Optional numeric vector of test survival times.
#' @param test_delta Optional numeric vector of test event indicators.
#' @param test_stratum Optional vector of test stratum membership.
#' @param criteria Character string: \code{"loss"} or \code{"CIndex"}.
#' @param ... Additional arguments (ignored).
#' 
#' @details
#' When \code{criteria = "loss"} and no test data are supplied, the plotted values are
#' \code{-2 * x$likelihood} (no normalization). When test data are provided,
#' performance is computed via \code{test_eval(..., criteria)}. The x-axis is shown
#' in decreasing \code{lambda} with a reversed log10 scale.
#' 
#' @return A \code{ggplot} object showing the performance curve.
#' @examples
#' data(ExampleData_highdim) 
#' 
#' train_dat_highdim <- ExampleData_highdim$train
#' test_dat_highdim <- ExampleData_highdim$test
#' beta_external_highdim <- ExampleData_highdim$beta_external
#' 
#' 
#' model_ridge <- coxkl_ridge(z = train_dat_highdim$z,
#'                            delta = train_dat_highdim$status,
#'                            time = train_dat_highdim$time,
#'                            beta = beta_external_highdim,
#'                            eta = 1)
#'
#' plot(
#'   model_ridge,
#'   test_z       = test_dat_highdim$z,
#'   test_time    = test_dat_highdim$time,
#'   test_delta   = test_dat_highdim$status,
#'   test_stratum = test_dat_highdim$stratum,
#'   criteria     = "CIndex"
#' )
#' 
#' @importFrom ggplot2 ggplot geom_segment aes geom_line geom_point labs coord_cartesian scale_x_reverse theme_minimal theme element_blank element_line element_text
#' @importFrom grid unit
#' @importFrom rlang .data
#' @exportS3Method plot coxkl_ridge
plot.coxkl_ridge <- function(x, test_z = NULL, test_time = NULL, test_delta = NULL,
                             test_stratum = NULL, criteria = c("loss", "CIndex"), ...) {
  criteria <- match.arg(criteria)
  if (!inherits(x, "coxkl_ridge")) stop("'x' must be of class 'coxkl_ridge'.", call. = FALSE)
  
  lambdas       <- x$lambda
  beta_mat      <- x$beta
  z_train       <- x$data$z
  time_train    <- x$data$time
  delta_train   <- x$data$delta
  stratum_train <- x$data$stratum
  
  using_train <- is.null(test_z) && is.null(test_time) && is.null(test_delta) && is.null(test_stratum)
  if (using_train) {
    test_z       <- z_train
    test_time    <- time_train
    test_delta   <- delta_train
    test_stratum <- stratum_train
  }
  
  n_eval <- nrow(as.matrix(test_z))
  
  if (criteria == "loss") {
    if (using_train) {
      metrics <- (-2 * x$likelihood) / n_eval
    } else {
      raw_metrics <- sapply(seq_along(lambdas), function(i)
        test_eval(
          test_z       = test_z,
          test_delta   = test_delta,
          test_time    = test_time,
          test_stratum = test_stratum,
          betahat      = beta_mat[, i],
          criteria     = "loss"
        )
      )
      metrics <- as.numeric(raw_metrics) / n_eval
    }
    opt_idx <- which.min(metrics)
  } else {
    metrics <- sapply(seq_along(lambdas), function(i)
      test_eval(
        test_z       = test_z,
        test_delta   = test_delta,
        test_time    = test_time,
        test_stratum = test_stratum,
        betahat      = beta_mat[, i],
        criteria     = "CIndex"
      )
    )
    opt_idx <- which.max(metrics)
  }
  
  df <- data.frame(lambda = lambdas, metric = as.numeric(metrics))
  df <- df[order(df$lambda, decreasing = TRUE), ]
  
  opt_lambda <- df$lambda[opt_idx]
  ylow  <- min(df$metric, na.rm = TRUE) * 0.995
  yhigh <- max(df$metric, na.rm = TRUE) * 1.005
  
  ggplot(df, aes(x = .data$lambda, y = .data$metric)) +
    geom_line(size = 1, color = "#7570B3") +
    geom_point(size = 2, color = "#7570B3") +
    geom_segment(
      data = data.frame(x = opt_lambda),
      aes(x = .data$x, xend = .data$x, y = ylow, yend = yhigh),
      inherit.aes = FALSE, color = "#D95F02", linewidth = 1, linetype = "dashed"
    ) +
    scale_x_reverse(trans = "log10") +
    labs(x = expression(lambda), y = if (criteria == "CIndex") "C Index" else "Loss") +
    coord_cartesian(ylim = c(ylow, yhigh)) +
    theme_biometrics()
}






#' Plot Model Performance vs Lambda for `coxkl_enet`
#'
#' @description
#' Plots model performance across the \code{lambda} sequence. Performance is
#' loss (\code{-2} times partial log-likelihood) or concordance index (C-index).
#' If no test data are provided, the curve uses the training data stored in \code{x$data}.
#'
#' @param x A fitted model object of class \code{"coxkl_enet"}.
#' @param test_z Optional numeric matrix of test covariates.
#' @param test_time Optional numeric vector of test survival times.
#' @param test_delta Optional numeric vector of test event indicators.
#' @param test_stratum Optional vector of test stratum membership.
#' @param criteria Character string: \code{"loss"} or \code{"CIndex"}.
#' @param ... Additional arguments (ignored).
#' 
#' 
#' @details
#' When \code{criteria = "loss"} and no test data are supplied, the plotted values are
#' \code{-2 * x$likelihood} (no normalization). When test data are provided,
#' performance is computed via \code{test_eval(..., criteria)}. The x-axis is shown
#' in decreasing \code{lambda} with a reversed log10 scale.
#'
#' @return A \code{ggplot} object showing the performance curve.
#' 
#' @examples
#' data(ExampleData_highdim) 
#' 
#' train_dat_highdim <- ExampleData_highdim$train
#' test_dat_highdim <- ExampleData_highdim$test
#' beta_external_highdim <- ExampleData_highdim$beta_external
#' 
#' model_enet <- coxkl_enet(z = train_dat_highdim$z,
#'                          delta = train_dat_highdim$status,
#'                          time = train_dat_highdim$time,
#'                          beta = beta_external_highdim,
#'                          eta = 1,
#'                          alpha = 1.0)
#' plot(model_enet,
#'      test_z = test_dat_highdim$z,
#'      test_time = test_dat_highdim$time,
#'      test_delta = test_dat_highdim$status,
#'      test_stratum = test_dat_highdim$stratum,
#'      criteria = "loss")
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_line geom_point labs coord_cartesian scale_x_reverse theme_minimal theme element_blank element_line element_text
#' @importFrom grid unit
#' @exportS3Method plot coxkl_enet
plot.coxkl_enet <- function(x, test_z = NULL, test_time = NULL, test_delta = NULL,
                            test_stratum = NULL, criteria = c("loss", "CIndex"), ...) {
  criteria <- match.arg(criteria)
  if (!inherits(x, "coxkl_enet")) stop("'x' must be of class 'coxkl_enet'.", call. = FALSE)
  
  lambdas       <- x$lambda
  beta_mat      <- x$beta
  z_train       <- x$data$z
  time_train    <- x$data$time
  delta_train   <- x$data$delta
  stratum_train <- x$data$stratum
  
  using_train <- is.null(test_z) && is.null(test_time) && is.null(test_delta) && is.null(test_stratum)
  if (using_train) {
    test_z       <- z_train
    test_time    <- time_train
    test_delta   <- delta_train
    test_stratum <- stratum_train
  }
  
  n_eval <- nrow(as.matrix(test_z))
  
  if (criteria == "loss") {
    if (using_train) {
      metrics <- (-2 * x$likelihood) / n_eval
    } else {
      raw_metrics <- sapply(seq_along(lambdas), function(i)
        test_eval(
          test_z       = test_z,
          test_delta   = test_delta,
          test_time    = test_time,
          test_stratum = test_stratum,
          betahat      = beta_mat[, i],
          criteria     = "loss"
        )
      )
      metrics <- as.numeric(raw_metrics) / n_eval
    }
    opt_idx <- which.min(metrics)
  } else {
    metrics <- sapply(seq_along(lambdas), function(i)
      test_eval(
        test_z       = test_z,
        test_delta   = test_delta,
        test_time    = test_time,
        test_stratum = test_stratum,
        betahat      = beta_mat[, i],
        criteria     = "CIndex"
      )
    )
    opt_idx <- which.max(metrics)
  }
  
  df <- data.frame(lambda = lambdas, metric = as.numeric(metrics))
  df <- df[order(df$lambda, decreasing = TRUE), ]
  
  opt_lambda <- df$lambda[opt_idx]
  ylow  <- min(df$metric, na.rm = TRUE) * 0.995
  yhigh <- max(df$metric, na.rm = TRUE) * 1.005
  
  ggplot(df, aes(x = .data$lambda, y = .data$metric)) +
    geom_line(size = 1, color = "#7570B3") +
    geom_point(size = 2, color = "#7570B3") +
    geom_segment(
      data = data.frame(x = opt_lambda),
      aes(x = .data$x, xend = .data$x, y = ylow, yend = yhigh),
      inherit.aes = FALSE, color = "#D95F02", linewidth = 1, linetype = "dashed"
    ) +
    scale_x_reverse(trans = "log10") +
    labs(x = expression(lambda), y = if (criteria == "CIndex") "C Index" else "Loss") +
    coord_cartesian(ylim = c(ylow, yhigh)) +
    theme_biometrics()
}


#' Internal theme for survkl plotting
#'
#' @keywords internal
#' @noRd
theme_biometrics <- function() {
  theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.length = unit(0.1, "cm"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 14, hjust = 0.0),
      legend.position = "none"
    )
}

