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
  
  etas <- object$eta
  beta_mat <- object$beta
  z_train <- object$data$z
  time_train <- object$data$time
  delta_train <- object$data$delta
  stratum_train <- object$data$stratum
  RS_external <- object$data$RS
  
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
  
  ext_perf <- NA_real_
  if (!is.null(RS_external)) {
    ext_perf <- test_eval(
      test_RS      = RS_external,
      test_delta   = test_delta,
      test_time    = test_time,
      test_stratum = test_stratum,
      criteria     = criteria
    )
  }
  
  df <- data.frame(eta = etas, metric = metrics)
  ylab <- if (criteria == "CIndex") "C Index" else "Loss"
  
  p <- ggplot(df, aes(x = eta, y = metric)) +
    geom_line(size = 1, color = "#7570B3") +
    geom_point(size = 2, color = "#7570B3") +
    labs(x = expression(eta), y = ylab) +
    theme_biometrics()
  
  if (is.finite(ext_perf)) {
    p <- p +
      geom_segment(
        inherit.aes = FALSE,
        x = 0, xend = max(df$eta, na.rm = TRUE),
        y = ext_perf, yend = ext_perf,
        color = "#D95F02", linetype = "dotted", linewidth = 1
      ) +
      coord_cartesian(
        ylim = c(min(c(df$metric, ext_perf), na.rm = TRUE) * 0.995,
                 max(c(df$metric, ext_perf), na.rm = TRUE) * 1.005)
      )
  } else {
    p <- p +
      coord_cartesian(
        ylim = c(min(df$metric, na.rm = TRUE) * 0.995,
                 max(df$metric, na.rm = TRUE) * 1.005)
      )
  }
  
  p
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
    scale_x_reverse(trans = "log10") +  # ✅ 单一 log10+reverse
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













