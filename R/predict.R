#' Predict Linear Predictors from a `coxkl` Object
#'
#' @description
#' Computes linear predictors for new data based on a fitted \code{coxkl} model.
#' If \code{eta} is supplied, predictions are returned for those \code{eta} values;
#' otherwise predictions are returned for all fitted \code{eta}s. Linear interpolation 
#' is applied if an intermediate \code{eta} value is requested.
#'
#' @param object A fitted model object of class \code{"coxkl"}.
#' @param newz A numeric matrix or data frame of new covariates (must match the
#'   dimension of the training design matrix used to fit the model).
#' @param eta Optional numeric vector of \code{eta} value(s) for which to predict.
#'   If \code{NULL}, predictions for all fitted \code{eta} values are returned.
#' @param ... Additional arguments.
#' 
#' @details
#' The linear predictors are computed as \code{as.matrix(newz) \%*\% beta}.
#' 
#' @return
#' A numeric matrix of linear predictors with one column per \code{eta} (sorted ascending).
#'
#' @seealso \code{\link{coef.coxkl}}
#'
#' @importFrom stats coef
#' @exportS3Method predict coxkl
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
#' Computes linear predictors for new data using a fitted \code{coxkl_ridge} model.
#' If \code{lambda} is supplied, predictions are returned for those \code{lambda}
#' values; otherwise predictions are returned for all fitted \code{lambda}s. When a
#' requested \code{lambda} lies between fitted values, coefficients are linearly
#' interpolated.
#'
#' @param object A fitted model object of class \code{"coxkl_ridge"}.
#' @param newz A numeric matrix or data frame of new covariates (same columns as in training data).
#' @param lambda Optional numeric value(s) specifying the regularization parameter(s)
#'   for which to predict. If `NULL`, predictions for all fitted \code{lambda} values are returned.
#' @param ... Additional arguments.
#' 
#' @details
#' The linear predictors are computed as \code{as.matrix(newz) \%*\% beta}.
#' 
#' @return
#' A numeric matrix of linear predictors.  
#' Each column corresponds to one \code{lambda}, sorted in descending order.
#'
#' @seealso \code{\link{coef.coxkl_ridge}}
#'
#' @importFrom stats coef
#' @exportS3Method predict coxkl_ridge
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
#' Computes linear predictors for new data using a fitted \code{coxkl_enet} model.
#' If \code{lambda} is supplied, predictions are returned for those \code{lambda}
#' values; otherwise predictions are returned for all fitted \code{lambda}s. When a
#' requested \code{lambda} lies between fitted values, coefficients are linearly
#' interpolated.
#'
#' @param object A fitted model object of class \code{"coxkl_enet"}.
#' @param newz A numeric matrix or data frame of new covariates (same columns as in training data).
#' @param lambda Optional numeric value(s) specifying the regularization parameter(s)
#'   for which to predict. If `NULL`, predictions for all fitted \code{lambda} values are returned.
#' @param ... Additional arguments.
#'
#' @details
#' The linear predictors are computed as \code{as.matrix(newz) \%*\% beta}.
#' 
#' @return
#' A numeric matrix of linear predictors.  
#' Each column corresponds to one \code{lambda}, sorted in descending order.
#'
#' @seealso \code{\link{coef.coxkl_enet}}
#' @importFrom stats coef
#' @exportS3Method predict coxkl_enet
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