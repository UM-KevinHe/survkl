#' Cox Proportional Hazards Model with KL Divergence for Data Integration
#' 
#' Fits a Cox proportional hazards model that incorporates external information
#' via a Kullback–Leibler (KL) divergence penalty. External information can be
#' supplied either as external risk scores (`RS`) or as external coefficients
#' (`beta`). The tuning parameter(s) `etas` control the strength of integration.
#'
#' @param z Numeric matrix of covariates with rows representing observations and
#'   columns representing predictor variables. All covariates must be numeric.
#' @param delta Numeric vector of event indicators (1 = event, 0 = censored).
#' @param time Numeric vector of observed event or censoring times. No sorting
#'   required.
#' @param stratum Optional numeric or factor vector defining strata.
#' @param RS Optional numeric vector or matrix of external risk scores. Length
#'   (or number of rows) must equal the number of observations. If not supplied,
#'   `beta` must be provided.
#' @param beta Optional numeric vector of external coefficients (e.g., from prior
#'   studies). Length must equal the number of columns in `z`. Use zeros to
#'   represent covariates without external information. If not supplied, `RS`
#'   must be provided.
#' @param etas Numeric vector of tuning parameters controlling the reliance on
#'   external information. Larger values place more weight on the external
#'   source.
#' @param tol Convergence tolerance for the optimization algorithm. Default is
#'   `1e-4`.
#' @param Mstop Maximum number of iterations for the optimization algorithm.
#'   Default is `100`.
#' @param backtrack Logical; if `TRUE`, backtracking line search is applied during
#'   optimization. Default is `FALSE`.
#' @param message Logical; if `TRUE`, progress messages are printed during model
#'   fitting. Default is `FALSE`.
#' @param data_sorted Logical; if `TRUE`, input data are assumed to be already
#'   sorted by stratum and time. Default is `FALSE`.
#' @param beta_initial Optional numeric vector of length `p` giving the starting
#'   value for the first `eta`. If `NULL`, a zero vector is used.
#'
#' @return
#' An object of class \code{"coxkl"} containing:
#' \itemize{
#'   \item \code{eta}: the fitted \eqn{\eta} sequence.
#'   \item \code{beta}: estimated coefficient matrix (\eqn{p \times |\eta|}).
#'   \item \code{linear.predictors}: matrix of linear predictors.
#'   \item \code{likelihood}: vector of partial likelihoods.
#'   \item \code{data}: a list containing the input data used in fitting
#'         (\code{z}, \code{time}, \code{delta}, \code{stratum}, \code{data_sorted}).
#' }
#' 
#' @details
#' If `beta` is supplied (length `ncol(z)`), external risk scores are computed
#' internally as `RS = z %*% beta`. If `RS` is supplied, it is used directly.
#' Data are optionally sorted by `stratum` (or a single stratum if `NULL`) and
#' increasing `time` when `data_sorted = FALSE`. Estimation proceeds over the
#' sorted data, and the returned `linear.predictors` are mapped back to the
#' original order. Optimization uses warm starts across the (ascending) `etas`
#' grid and supports backtracking line search when `backtrack = TRUE`.
#'
#' Internally, the routine computes a stratum-wise adjusted event indicator
#' (`delta_tilde`) and maximizes a KL-regularized partial likelihood. The current
#' implementation fixes `lambda = 0` in the low-level optimizer and exposes
#' `etas` as the primary tuning control.
#'
#' @examples
#' data(Exampledata_lowdim)
#' 
#' train_dat_lowdim <- ExampleData_lowdim$train
#' beta_external_good_lowdim <- ExampleData_lowdim$beta_external_good
#' 
#' model <- coxkl(z = train_dat_lowdim$z,
#'      delta = train_dat_lowdim$status,
#'      time = train_dat_lowdim$time,
#'      stratum = train_dat_lowdim$stratum,
#'      RS = NULL,
#'      beta = beta_external_good_lowdim,
#'      etas = c(0:5))
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @useDynLib survkl, .registration = TRUE
#'
#' @export

coxkl <- function(z, delta, time, stratum = NULL,
                  RS = NULL, beta = NULL, 
                  etas, tol = 1.0e-4, Mstop = 100,
                  backtrack = FALSE,
                  message = FALSE,
                  data_sorted = FALSE,
                  beta_initial = NULL){
  
  if (is.null(RS) && is.null(beta)) {
    stop("No external information is provided. Either RS or beta must be provided.")
  } else if (is.null(RS) && !is.null(beta)) {
    if (length(beta) == ncol(z)) {
      if (message) message("External beta information is used.")
      RS <- as.matrix(z) %*% as.matrix(beta)
    } else {
      stop("The dimension of beta does not match the number of columns in z.")
    }
  } else if (!is.null(RS)) {
    RS <- as.matrix(RS)
    if (message) message("External Risk Score information is used.")
  }
  
  input_data <- list(z = z, time = time, delta = delta, stratum = stratum, RS = RS)
  
  if (!data_sorted) {
    ## ---- Sorting Section ----
    if (is.null(stratum)) {
      if (message) warning("Stratum information not provided. All data is assumed to originate from a single stratum!", call. = FALSE)
      stratum <- rep(1, nrow(z))
    } else {
      stratum <- match(stratum, unique(stratum))
    }
    time_order <- order(stratum, time)
    time <- as.numeric(time[time_order])
    stratum <- as.numeric(stratum[time_order])
    z_mat <- as.matrix(z)[time_order, , drop = FALSE]
    delta <- as.numeric(delta[time_order])
    RS <- as.numeric(RS[time_order, , drop = FALSE])
  } else {
    z_mat <- as.matrix(z)
    time <- as.numeric(time)
    delta <- as.numeric(delta)
    stratum <- as.numeric(stratum)
    RS <- as.numeric(RS)
  }
  
  etas <- sort(etas)
  n_eta <- length(etas)
  LP_mat <- matrix(NA, nrow = nrow(z_mat), ncol = n_eta)
  beta_mat <- matrix(NA, nrow = ncol(z_mat), ncol = n_eta)
  likelihood_mat <- rep(NA, n_eta)
  
  eta_names <- round(etas, 4)
  colnames(LP_mat) <- eta_names
  colnames(beta_mat) <- eta_names
  names(likelihood_mat) <- eta_names
  
  n.each_stratum <- as.numeric(table(stratum))
  delta_tilde <- calculateDeltaTilde(delta, time, RS, n.each_stratum)
  
  if (is.null(beta_initial)){
    beta_initial <- rep(0, ncol(z_mat))
  }
  
  if (message) {
    cat("Cross-validation over eta sequence:\n")
    pb <- txtProgressBar(min = 0, max = n_eta, style = 3, width = 30)
  }
  
  for (i in seq_along(etas)){  #"etas" already in ascending order
    eta <- etas[i]
    beta_est <- KL_Cox_Estimate_cpp(z_mat, delta, delta_tilde, n.each_stratum, eta, beta_initial,
                                    tol, Mstop, lambda = 0, backtrack = backtrack, message = F)
    LP <- z_mat %*% as.matrix(beta_est)
    LP_mat[, i] <- LP
    beta_mat[, i] <- beta_est
    likelihood_mat[i] <- pl_cal_theta(LP, delta, n.each_stratum)
    
    beta_initial <- beta_est  # "warm start"
    if (message) setTxtProgressBar(pb, i)
  }
  if (message) close(pb)
  
  if (data_sorted == FALSE){
    LinPred_original <- matrix(NA_real_, nrow = length(time_order), ncol = n_eta)
    LinPred_original[time_order, ] <- LP_mat
  } else {
    LinPred_original <- LP_mat
  }
  
  structure(list(
    eta = etas,
    beta = beta_mat,
    linear.predictors = LinPred_original,
    likelihood = likelihood_mat,
    data = input_data
  ), class = "coxkl")
}

