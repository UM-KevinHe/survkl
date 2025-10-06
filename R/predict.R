#' Predict Linear Predictors from a coxkl Object
#'
#' @description
#' Computes linear predictors for new data based on a fitted \code{coxkl} model.
#' Users can specify one or more \code{eta} values; if not provided, predictions
#' are returned for all fitted \code{eta} values. Linear interpolation is applied
#' if an intermediate \code{eta} value is requested.
#'
#' @param object A fitted model object of class \code{"coxkl"}.
#' @param newz A numeric matrix or data frame of new covariates (must match the
#'   dimension of the training design matrix used to fit the model).
#' @param eta Optional numeric value(s) specifying which \code{eta} values to use
#'   for prediction. If `NULL`, predictions for all fitted \code{eta} values are returned.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A numeric matrix of linear predictors.  
#' Each column corresponds to one value of \code{eta}, sorted in ascending order.
#'
#' @export
predict.coxkl <- function(object, newz, eta = NULL, ...) {
  if (!inherits(object, "coxkl")) stop("'object' must be of class 'coxkl'.", call. = FALSE)
  if (missing(newz)) stop("Argument 'newz' must be provided.", call. = FALSE)
  
  beta_mat <- coef(object, eta = eta)  # handles sorting & interpolation
  Z <- as.matrix(newz)
  
  if (nrow(Z) == 0L || ncol(Z) != nrow(beta_mat)) {
    stop("Dimension mismatch between 'newz' and model coefficients.", call. = FALSE)
  }
  
  lp <- Z %*% beta_mat
  colnames(lp) <- colnames(beta_mat)
  return(lp)
}


#' Predict Linear Predictors from a coxkl_ridge Object
#'
#' @description
#' Computes linear predictors for new data from a fitted \code{coxkl_ridge} model.
#' Users can specify one or more \code{lambda} values; if not provided, predictions
#' are returned for all fitted \code{lambda} values. Linear interpolation is applied
#' if an intermediate \code{lambda} value is requested.
#'
#' @param object A fitted model object of class \code{"coxkl_ridge"}.
#' @param newz A numeric matrix or data frame of new covariates (same columns as in training data).
#' @param lambda Optional numeric value(s) specifying the regularization parameter(s)
#'   for which to predict. If `NULL`, predictions for all fitted \code{lambda} values are returned.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A numeric matrix of linear predictors.  
#' Each column corresponds to one \code{lambda}, sorted in descending order.
#'
#' @export
predict.coxkl_ridge <- function(object, newz, lambda = NULL, ...) {
  if (!inherits(object, "coxkl_ridge")) stop("'object' must be of class 'coxkl_ridge'.", call. = FALSE)
  if (missing(newz)) stop("Argument 'newz' must be provided.", call. = FALSE)
  
  beta_mat <- coef(object, lambda = lambda)  # handles descending sort & interpolation
  Z <- as.matrix(newz)
  
  if (nrow(Z) == 0L || ncol(Z) != nrow(beta_mat)) {
    stop("Dimension mismatch between 'newz' and model coefficients.", call. = FALSE)
  }
  
  lp <- Z %*% beta_mat
  colnames(lp) <- colnames(beta_mat)
  return(lp)
}



#' Predict Linear Predictors from a coxkl_enet Object
#'
#' @description
#' Computes linear predictors for new data from a fitted \code{coxkl_enet} model.
#' Users can specify one or more \code{lambda} values; if not provided, predictions
#' are returned for all fitted \code{lambda} values. Linear interpolation is applied
#' if an intermediate \code{lambda} value is requested.
#'
#' @param object A fitted model object of class \code{"coxkl_enet"}.
#' @param newz A numeric matrix or data frame of new covariates (same columns as in training data).
#' @param lambda Optional numeric value(s) specifying the regularization parameter(s)
#'   for which to predict. If `NULL`, predictions for all fitted \code{lambda} values are returned.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A numeric matrix of linear predictors.  
#' Each column corresponds to one \code{lambda}, sorted in descending order.
#'
#' @export
predict.coxkl_enet <- function(object, newz, lambda = NULL, ...) {
  if (!inherits(object, "coxkl_enet")) stop("'object' must be of class 'coxkl_enet'.", call. = FALSE)
  if (missing(newz)) stop("Argument 'newz' must be provided.", call. = FALSE)
  
  beta_mat <- coef(object, lambda = lambda)  # handles descending sort & interpolation
  Z <- as.matrix(newz)
  
  if (nrow(Z) == 0L || ncol(Z) != nrow(beta_mat)) {
    stop("Dimension mismatch between 'newz' and model coefficients.", call. = FALSE)
  }
  
  lp <- Z %*% beta_mat
  colnames(lp) <- colnames(beta_mat)
  return(lp)
}

