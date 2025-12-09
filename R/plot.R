#' Plot Model Performance vs Eta for `coxkl`
#'
#' @description
#' Plots model performance across the \code{eta} sequence. Performance is either
#' loss (\code{-2} times partial log-likelihood) or concordance index (C-index).
#' If no test data are provided, the curve is computed on the training data stored
#' in \code{object$data}.
#'
#' @param object A fitted model object of class \code{"coxkl"}.
#' @param test_z Optional numeric matrix of test covariates.
#' @param test_time Optional numeric vector of test survival times.
#' @param test_delta Optional numeric vector of test event indicators.
#' @param test_stratum Optional vector of test stratum membership.
#' @param criteria Character string: \code{"loss"} or \code{"CIndex"}.
#' @param ... Additional arguments (ignored).
#' 
#' @details
#' When \code{criteria = "loss"} and no test data are supplied, the plotted values are
#' \code{(-2 * object$likelihood) / n}, where \code{n} is the number of rows in the
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
#' data(Exampledata_lowdim)
#' 
#' train_dat_lowdim <- ExampleData_lowdim$train
#' beta_external_good_lowdim <- ExampleData_lowdim$beta_external_good
#' 
#' model <- coxkl(z = train_dat_lowdim$z,
#'                delta = train_dat_lowdim$status,
#'                time = train_dat_lowdim$time,
#'                stratum = train_dat_lowdim$stratum,
#'                RS = NULL,
#'                beta = beta_external_good_lowdim,
#'                etas = c(0:5))
#'      
#' plot(model)      
#'      
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_segment labs coord_cartesian
#' @importFrom ggplot2 theme_minimal theme element_blank element_line element_text
#' @importFrom grid unit
#' @exportS3Method plot coxkl
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
  
  n_eval <- nrow(as.matrix(test_z))
  
  if (criteria == "loss") {
    if (using_train) {
      metrics <- (-2 * object$likelihood) / n_eval
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




#' Plot Model Performance vs Lambda for `coxkl_ridge`
#'
#' @description
#' Plots model performance across the \code{lambda} sequence. Performance is
#' loss (\code{-2} times partial log-likelihood) or concordance index (C-index).
#' If no test data are provided, the curve uses the training data stored in \code{object$data}.
#'
#' @param object A fitted model object of class \code{"coxkl_ridge"}.
#' @param test_z Optional numeric matrix of test covariates.
#' @param test_time Optional numeric vector of test survival times.
#' @param test_delta Optional numeric vector of test event indicators.
#' @param test_stratum Optional vector of test stratum membership.
#' @param criteria Character string: \code{"loss"} or \code{"CIndex"}.
#' @param ... Additional arguments (ignored).
#' 
#' @details
#' When \code{criteria = "loss"} and no test data are supplied, the plotted values are
#' \code{-2 * object$likelihood} (no normalization). When test data are provided,
#' performance is computed via \code{test_eval(..., criteria)}. The x-axis is shown
#' in decreasing \code{lambda} with a reversed log10 scale.
#' 
#' data(ExampleData_highdim) 
#' 
#' train_dat_highdim <- ExampleData_highdim$train
#' beta_external_highdim <- ExampleData_highdim$beta_external
#' 
#' model_ridge <- coxkl_ridge(z = train_dat_highdim$z,
#'                            delta = train_dat_highdim$status,
#'                            time = train_dat_highdim$time,
#'                            stratum = NULL,
#'                            RS = NULL,
#'                            beta = beta_external_highdim,
#'                            message = TRUE)
#'  
#' plot(model_ridge)
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs coord_cartesian scale_x_reverse
#' @importFrom ggplot2 theme_minimal theme element_blank element_line element_text
#' @importFrom grid unit
#' @exportS3Method plot coxkl_ridge
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




#' Plot Model Performance vs Lambda for `coxkl_enet`
#'
#' @description
#' Plots model performance across the \code{lambda} sequence. Performance is
#' loss (\code{-2} times partial log-likelihood) or concordance index (C-index).
#' If no test data are provided, the curve uses the training data stored in \code{object$data}.
#'
#' @param object A fitted model object of class \code{"coxkl_enet"}.
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
#' \code{-2 * object$likelihood} (no normalization). When test data are provided,
#' performance is computed via \code{test_eval(..., criteria)}. The x-axis is shown
#' in decreasing \code{lambda} with a reversed log10 scale.
#'
#' @return A \code{ggplot} object showing the performance curve.
#' 
#' @examples
#' data(ExampleData_highdim) 
#' 
#' train_dat_highdim <- ExampleData_highdim$train
#' beta_external_highdim <- ExampleData_highdim$beta_external
#' 
#' model_enet <- coxkl_enet(z = train_dat_highdim$z,
#'                          delta = train_dat_highdim$status,
#'                          time = train_dat_highdim$time,
#'                          stratum = NULL,
#'                          RS = NULL,
#'                          beta = beta_external_highdim,
#'                          eta = 0,
#'                          alpha = 1.0,
#'                          message = TRUE)
#'                          
#' plot(model_enet)
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs coord_cartesian scale_x_reverse
#' @importFrom ggplot2 theme_minimal theme element_blank element_line element_text
#' @importFrom grid unit
#' @exportS3Method plot coxkl_enet
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


#' Theme helper for biometrics-style plots (internal)
#'
#' @description
#' Internal \code{ggplot2} theme used by plotting methods in this package.
#'
#' @return A \code{ggplot2} theme object.
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






