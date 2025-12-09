#' Extract Coefficients from a `coxkl` Object
#' 
#' @description
#' Extracts the estimated regression coefficients (\code{beta}) from a fitted
#' \code{coxkl} object. Optionally, a value (or vector) of \code{eta} can be
#' supplied. If the requested \code{eta} values are not in the fitted sequence,
#' linear interpolation is performed between the nearest neighboring \code{eta}
#' values; out-of-range requests error.
#'
#' @param object An object of class \code{"coxkl"}, typically the result of
#'   \code{\link{coxkl}}.
#' @param eta Optional numeric value or vector specifying the \eqn{\eta}
#'   values for which to extract (or interpolate) coefficients. If \code{NULL},
#'   all estimated coefficients are returned.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A numeric matrix of regression coefficients.  
#' Each column corresponds to one value of \code{eta}, sorted in ascending order.
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
#' coef(model)
#'
#' @exportS3Method coef coxkl
coef.coxkl <- function(object, eta = NULL, ...) {
  if (!inherits(object, "coxkl")) {
    stop("'object' must be of class 'coxkl'.", call. = FALSE)
  }
  
  etas <- object$eta
  beta_mat <- object$beta
  
  if (is.null(eta)) {
    # Return all betas sorted by eta
    ord <- order(etas)
    return(beta_mat[, ord, drop = FALSE])
  }
  
  if (!is.numeric(eta)) {
    stop("'eta' must be numeric.", call. = FALSE)
  }
  
  eta_min <- min(etas)
  eta_max <- max(etas)
  
  if (any(eta < eta_min | eta > eta_max)) {
    stop(paste0("Some 'eta' values are outside the fitted range [", 
                eta_min, ", ", eta_max, "]."), call. = FALSE)
  }
  
  # Sort eta before processing for consistent output order
  eta_sorted <- sort(unique(eta))
  beta_interp <- matrix(NA_real_, nrow = nrow(beta_mat), ncol = length(eta_sorted))
  colnames(beta_interp) <- as.character(eta_sorted)
  rownames(beta_interp) <- rownames(beta_mat)
  
  for (j in seq_along(eta_sorted)) {
    target_eta <- eta_sorted[j]
    
    if (target_eta %in% etas) {
      beta_interp[, j] <- beta_mat[, which(etas == target_eta)]
    } else {
      lower_idx <- max(which(etas < target_eta))
      upper_idx <- min(which(etas > target_eta))
      
      eta_lower <- etas[lower_idx]
      eta_upper <- etas[upper_idx]
      w <- (target_eta - eta_lower) / (eta_upper - eta_lower)
      
      beta_interp[, j] <- (1 - w) * beta_mat[, lower_idx] + w * beta_mat[, upper_idx]
      
      warning(sprintf(
        "Linear interpolation performed between eta = %.3f and eta = %.3f for eta = %.3f.",
        eta_lower, eta_upper, target_eta
      ), call. = FALSE)
    }
  }
  
  return(beta_interp)
}


#' Extract Coefficients from a `coxkl_ridge` Object
#' 
#' @description
#' Extracts the estimated regression coefficients (\code{beta}) from a fitted
#' \code{coxkl_ridge} object. Optionally, one or more \code{lambda} values can be
#' supplied. If requested \code{lambda} values are not in the fitted sequence,
#' linear interpolation is performed between nearest neighbors; out-of-range
#' requests error.
#'
#' @param object An object of class \code{"coxkl_ridge"}, typically the result of
#'   \code{\link{coxkl_ridge}}.
#' @param lambda Optional numeric value or vector specifying the regularization
#'   parameter(s) for which to extract (or interpolate) coefficients. If \code{NULL},
#'   all estimated coefficients are returned.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A numeric matrix of regression coefficients.  
#' Each column corresponds to one value of \code{lambda}, sorted in \emph{descending} order.
#' 
#' @examples
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
#'                            message = T)
#' coef(model_ridge)
#'
#' @exportS3Method coef coxkl_ridge
coef.coxkl_ridge <- function(object, lambda = NULL, ...) {
  if (!inherits(object, "coxkl_ridge")) {
    stop("'object' must be of class 'coxkl_ridge'.", call. = FALSE)
  }
  
  lambdas <- object$lambda
  beta_mat <- object$beta
  
  # If no lambda specified, return all in descending order
  if (is.null(lambda)) {
    ord <- order(lambdas, decreasing = TRUE)
    return(beta_mat[, ord, drop = FALSE])
  }
  
  if (!is.numeric(lambda)) {
    stop("'lambda' must be numeric.", call. = FALSE)
  }
  
  lambda_min <- min(lambdas)
  lambda_max <- max(lambdas)
  
  if (any(lambda < lambda_min | lambda > lambda_max)) {
    stop(paste0("Some 'lambda' values are outside the fitted range [", 
                lambda_min, ", ", lambda_max, "]."), call. = FALSE)
  }
  
  # Sort in descending order (for output)
  lambda_sorted <- sort(unique(lambda), decreasing = TRUE)
  beta_interp <- matrix(NA_real_, nrow = nrow(beta_mat), ncol = length(lambda_sorted))
  colnames(beta_interp) <- as.character(lambda_sorted)
  rownames(beta_interp) <- rownames(beta_mat)
  
  for (j in seq_along(lambda_sorted)) {
    target_lambda <- lambda_sorted[j]
    
    if (target_lambda %in% lambdas) {
      beta_interp[, j] <- beta_mat[, which(lambdas == target_lambda)]
    } else {
      lower_idx <- max(which(lambdas < target_lambda))
      upper_idx <- min(which(lambdas > target_lambda))
      
      lambda_lower <- lambdas[lower_idx]
      lambda_upper <- lambdas[upper_idx]
      w <- (target_lambda - lambda_lower) / (lambda_upper - lambda_lower)
      
      beta_interp[, j] <- (1 - w) * beta_mat[, lower_idx] + w * beta_mat[, upper_idx]
      
      warning(sprintf(
        "Linear interpolation performed between lambda = %.3e and lambda = %.3e for lambda = %.3e.",
        lambda_lower, lambda_upper, target_lambda
      ), call. = FALSE)
    }
  }
  
  return(beta_interp)
}



#' Extract Coefficients from a `coxkl_enet` Object
#'
#' @description
#' Extracts the estimated regression coefficients (\code{beta}) from a fitted
#' \code{coxkl_enet} object. Optionally, one or more \code{lambda} values can be
#' supplied. If requested \code{lambda} values are not in the fitted sequence,
#' linear interpolation is performed between nearest neighbors; out-of-range
#' requests error.
#'
#' @param object An object of class \code{"coxkl_enet"}, typically the result of
#'   \code{\link{coxkl_enet}}.
#' @param lambda Optional numeric value or vector specifying the regularization
#'   parameter(s) for which to extract (or interpolate) coefficients. If \code{NULL},
#'   all estimated coefficients are returned.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A numeric matrix of regression coefficients; each column corresponds to one
#' value of \code{lambda}, sorted in \emph{descending} order.
#' 
#' @examples
#' data(ExampleData_highdim) 
#' 
#' train_dat_highdim <- ExampleData_highdim$train
#' beta_external_highdim <- ExampleData_highdim$beta_external
#' 
#' enet_model <- coxkl_enet(z = train_dat_highdim$z,
#'                          delta = train_dat_highdim$status,
#'                          time = train_dat_highdim$time,
#'                          stratum = NULL,
#'                          RS = NULL,
#'                          beta = beta_external_highdim,
#'                          eta = 0,
#'                          alpha = 1.0,
#'                          message = T)
#' coef(enet_model)                         
#'    
#' @exportS3Method coef coxkl_enet
coef.coxkl_enet <- function(object, lambda = NULL, ...) {
  if (!inherits(object, "coxkl_enet")) {
    stop("'object' must be of class 'coxkl_enet'.", call. = FALSE)
  }
  
  lambdas <- object$lambda
  beta_mat <- object$beta
  
  if (is.null(lambda)) {
    ord <- order(lambdas, decreasing = TRUE)
    return(beta_mat[, ord, drop = FALSE])
  }
  
  if (!is.numeric(lambda)) {
    stop("'lambda' must be numeric.", call. = FALSE)
  }
  
  lambda_min <- min(lambdas)
  lambda_max <- max(lambdas)
  
  if (any(lambda < lambda_min | lambda > lambda_max)) {
    stop(paste0("Some 'lambda' values are outside the fitted range [", 
                lambda_min, ", ", lambda_max, "]."), call. = FALSE)
  }
  
  lambda_sorted <- sort(unique(lambda), decreasing = TRUE)
  beta_interp <- matrix(NA_real_, nrow = nrow(beta_mat), ncol = length(lambda_sorted))
  colnames(beta_interp) <- as.character(lambda_sorted)
  rownames(beta_interp) <- rownames(beta_mat)
  
  for (j in seq_along(lambda_sorted)) {
    target_lambda <- lambda_sorted[j]
    
    if (target_lambda %in% lambdas) {
      beta_interp[, j] <- beta_mat[, which(lambdas == target_lambda)]
    } else {
      lower_idx <- max(which(lambdas < target_lambda))
      upper_idx <- min(which(lambdas > target_lambda))
      
      lambda_lower <- lambdas[lower_idx]
      lambda_upper <- lambdas[upper_idx]
      w <- (target_lambda - lambda_lower) / (lambda_upper - lambda_lower)
      
      beta_interp[, j] <- (1 - w) * beta_mat[, lower_idx] + w * beta_mat[, upper_idx]
      
      warning(sprintf(
        "Linear interpolation performed between lambda = %.3e and lambda = %.3e for lambda = %.3e.",
        lambda_lower, lambda_upper, target_lambda
      ), call. = FALSE)
    }
  }
  
  return(beta_interp)
}
