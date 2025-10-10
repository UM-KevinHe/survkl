#' Plot Model Performance vs Eta for coxkl
#'
#' @description
#' Plots model performance (loss or concordance index) across the \code{eta} sequence.
#' If no test data are provided, the plot is based on training data likelihoods.
#'
#' @param object A fitted model object of class \code{"coxkl"}.
#' @param test_z Optional numeric matrix of test covariates.
#' @param test_time Optional numeric vector of test survival times.
#' @param test_delta Optional numeric vector of test event indicators.
#' @param test_stratum Optional vector of test stratum membership.
#' @param criteria Character string: \code{"loss"} or \code{"CIndex"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot} object showing the performance curve.
#' @export
plot.coxkl <- function(object, test_z = NULL, test_time = NULL, test_delta = NULL,
                       test_stratum = NULL, criteria = c("loss", "CIndex"), ...) {
  criteria <- match.arg(criteria)
  if (!inherits(object, "coxkl")) stop("'object' must be of class 'coxkl'.", call. = FALSE)
  
  etas          <- object$eta
  beta_mat      <- object$beta
  z_train       <- object$data$z
  time_train    <- object$data$time
  delta_train   <- object$data$delta
  stratum_train <- object$data$stratum
  
  using_train <- is.null(test_z) && is.null(test_time) && is.null(test_delta) && is.null(test_stratum)
  if (using_train) {
    test_z       <- z_train
    test_time    <- time_train
    test_delta   <- delta_train
    test_stratum <- stratum_train
  }
  
  if (using_train && criteria == "loss") {
    metrics <- -2 * object$likelihood
  } else {
    metrics <- sapply(seq_along(etas), function(i)
      test_eval(
        test_z       = test_z,
        test_delta   = test_delta,
        test_time    = test_time,
        test_stratum = test_stratum,
        betahat      = beta_mat[, i],
        criteria     = criteria
      )
    )
  }
  
  idx0 <- which.min(abs(etas - 0))
  x0   <- etas[idx0]
  y0   <- as.numeric(metrics[idx0])
  xmax <- max(etas, na.rm = TRUE)
  
  df   <- data.frame(eta = etas, metric = as.numeric(metrics))
  ylab <- if (criteria == "CIndex") "C Index" else "Loss"
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = eta, y = metric)) +
    ggplot2::geom_line(size = 1, color = "#7570B3") +
    ggplot2::geom_point(size = 2, color = "#7570B3") +
    ggplot2::geom_segment(
      data = data.frame(x0 = x0, y0 = y0, xmax = xmax),
      ggplot2::aes(x = x0, xend = xmax, y = y0, yend = y0),
      inherit.aes = FALSE, color = "#1B9E77", linetype = "dotted", linewidth = 1
    ) +
    ggplot2::geom_point(
      data = data.frame(x0 = x0, y0 = y0),
      ggplot2::aes(x = x0, y = y0),
      inherit.aes = FALSE, color = "#1B9E77", shape = 16, size = 3
    ) +
    ggplot2::labs(x = expression(eta), y = ylab) +
    theme_biometrics() +
    ggplot2::coord_cartesian(
      ylim = c(min(c(df$metric, y0), na.rm = TRUE) * 0.995,
               max(c(df$metric, y0), na.rm = TRUE) * 1.005)
    )
  
  return(p)
}




#' Plot Model Performance vs Lambda for coxkl_ridge
#'
#' @description
#' Plots model performance (loss or concordance index) across the \code{lambda} sequence.
#' If no test data are provided, the plot uses training data likelihoods.
#'
#' @param object A fitted model object of class \code{"coxkl_ridge"}.
#' @param test_z Optional numeric matrix of test covariates.
#' @param test_time Optional numeric vector of test survival times.
#' @param test_delta Optional numeric vector of test event indicators.
#' @param test_stratum Optional vector of test stratum membership.
#' @param criteria Character string: \code{"loss"} or \code{"CIndex"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot} object showing the performance curve.
#' @export
plot.coxkl_ridge <- function(object, test_z = NULL, test_time = NULL, test_delta = NULL,
                             test_stratum = NULL, criteria = c("loss", "CIndex"), ...) {
  criteria <- match.arg(criteria)
  if (!inherits(object, "coxkl_ridge")) stop("'object' must be of class 'coxkl_ridge'.", call. = FALSE)
  
  lambdas <- object$lambda
  beta_mat <- object$beta
  z_train <- object$data$z
  time_train <- object$data$time
  delta_train <- object$data$delta
  stratum_train <- object$data$stratum
  
  using_train <- is.null(test_z) && is.null(test_time) && is.null(test_delta) && is.null(test_stratum)
  if (using_train) {
    test_z <- z_train
    test_time <- time_train
    test_delta <- delta_train
    test_stratum <- stratum_train
  }
  
  if (using_train && criteria == "loss") {
    metrics <- -2 * object$likelihood
  } else {
    metrics <- sapply(seq_along(lambdas), function(i)
      test_eval(
        test_z       = test_z,
        test_delta   = test_delta,
        test_time    = test_time,
        test_stratum = test_stratum,
        betahat      = beta_mat[, i],
        criteria     = criteria
      )
    )
  }
  
  df <- data.frame(lambda = lambdas, metric = metrics)
  df <- df[order(df$lambda, decreasing = TRUE), ]
  ylab <- if (criteria == "CIndex") "C Index" else "Loss"
  
  ggplot(df, aes(x = lambda, y = metric)) +
    geom_line(size = 1, color = "#7570B3") +
    geom_point(size = 2, color = "#7570B3") +
    scale_x_reverse(trans = "log10") +
    labs(x = expression(lambda), y = ylab) +
    coord_cartesian(
      ylim = c(min(df$metric, na.rm = TRUE) * 0.995,
               max(df$metric, na.rm = TRUE) * 1.005)
    ) +
    theme_biometrics()
}




#' Plot Model Performance vs Lambda for coxkl_enet
#'
#' @description
#' Plots model performance (loss or concordance index) across the \code{lambda} sequence.
#' If no test data are provided, the plot uses training data likelihoods.
#'
#' @param object A fitted model object of class \code{"coxkl_enet"}.
#' @param test_z Optional numeric matrix of test covariates.
#' @param test_time Optional numeric vector of test survival times.
#' @param test_delta Optional numeric vector of test event indicators.
#' @param test_stratum Optional vector of test stratum membership.
#' @param criteria Character string: \code{"loss"} or \code{"CIndex"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot} object showing the performance curve.
#' @export
plot.coxkl_enet <- function(object, test_z = NULL, test_time = NULL, test_delta = NULL,
                            test_stratum = NULL, criteria = c("loss", "CIndex"), ...) {
  criteria <- match.arg(criteria)
  if (!inherits(object, "coxkl_enet")) stop("'object' must be of class 'coxkl_enet'.", call. = FALSE)
  
  lambdas <- object$lambda
  beta_mat <- object$beta
  z_train <- object$data$z
  time_train <- object$data$time
  delta_train <- object$data$delta
  stratum_train <- object$data$stratum
  
  using_train <- is.null(test_z) && is.null(test_time) && is.null(test_delta) && is.null(test_stratum)
  if (using_train) {
    test_z <- z_train
    test_time <- time_train
    test_delta <- delta_train
    test_stratum <- stratum_train
  }
  
  if (using_train && criteria == "loss") {
    metrics <- -2 * object$likelihood
  } else {
    metrics <- sapply(seq_along(lambdas), function(i)
      test_eval(
        test_z       = test_z,
        test_delta   = test_delta,
        test_time    = test_time,
        test_stratum = test_stratum,
        betahat      = beta_mat[, i],
        criteria     = criteria
      )
    )
  }
  
  df <- data.frame(lambda = lambdas, metric = metrics)
  df <- df[order(df$lambda, decreasing = TRUE), ]
  ylab <- if (criteria == "CIndex") "C Index" else "Loss"
  
  ggplot(df, aes(x = lambda, y = metric)) +
    geom_line(size = 1, color = "#7570B3") +
    geom_point(size = 2, color = "#7570B3") +
    scale_x_reverse(trans = "log10") +
    labs(x = expression(lambda), y = ylab) +
    coord_cartesian(
      ylim = c(min(df$metric, na.rm = TRUE) * 0.995,
               max(df$metric, na.rm = TRUE) * 1.005)
    ) +
    theme_biometrics()
}






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













