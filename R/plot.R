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
  
  if (is.null(test_z)) {
    if (criteria == "loss") {
      metrics <- -2 * object$likelihood
    } else {
      test_z <- object$data$z
      test_time <- object$data$time
      test_delta <- object$data$delta
      test_stratum <- object$data$stratum
      metrics <- sapply(seq_along(etas), function(i)
        test_eval(test_z, test_delta, test_time, test_stratum,
                  betahat = beta_mat[, i], criteria = "CIndex"))
    }
  } else {
    metrics <- sapply(seq_along(etas), function(i)
      test_eval(test_z, test_delta, test_time, test_stratum,
                betahat = beta_mat[, i], criteria = criteria))
  }
  
  df <- data.frame(eta = etas, metric = metrics)
  ylab <- if (criteria == "CIndex") "Concordance Index" else "Loss"
  
  ggplot2::ggplot(df, ggplot2::aes(x = eta, y = metric)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = expression(eta), y = ylab)
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
  
  if (is.null(test_z)) {
    if (criteria == "loss") {
      metrics <- -2 * object$likelihood
    } else {
      test_z <- object$data$z
      test_time <- object$data$time
      test_delta <- object$data$delta
      test_stratum <- object$data$stratum
      metrics <- sapply(seq_along(lambdas), function(i)
        test_eval(test_z, test_delta, test_time, test_stratum,
                  betahat = beta_mat[, i], criteria = "CIndex"))
    }
  } else {
    metrics <- sapply(seq_along(lambdas), function(i)
      test_eval(test_z, test_delta, test_time, test_stratum,
                betahat = beta_mat[, i], criteria = criteria))
  }
  
  df <- data.frame(lambda = lambdas, metric = metrics)
  df <- df[order(df$lambda, decreasing = TRUE), ]
  ylab <- if (criteria == "CIndex") "Concordance Index" else "Loss"
  
  ggplot2::ggplot(df, ggplot2::aes(x = lambda, y = metric)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_x_reverse() +  # 让 log 轴从大→小
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = expression(lambda), y = ylab)
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
  
  if (is.null(test_z)) {
    if (criteria == "loss") {
      metrics <- -2 * object$likelihood
    } else {
      test_z <- object$data$z
      test_time <- object$data$time
      test_delta <- object$data$delta
      test_stratum <- object$data$stratum
      metrics <- sapply(seq_along(lambdas), function(i)
        test_eval(test_z, test_delta, test_time, test_stratum,
                  betahat = beta_mat[, i], criteria = "CIndex"))
    }
  } else {
    metrics <- sapply(seq_along(lambdas), function(i)
      test_eval(test_z, test_delta, test_time, test_stratum,
                betahat = beta_mat[, i], criteria = criteria))
  }
  
  df <- data.frame(lambda = lambdas, metric = metrics)
  df <- df[order(df$lambda, decreasing = TRUE), ]
  ylab <- if (criteria == "CIndex") "Concordance Index" else "Loss"
  
  ggplot2::ggplot(df, ggplot2::aes(x = lambda, y = metric)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_x_reverse() +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = expression(lambda), y = ylab)
}
