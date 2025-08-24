#' Calculate Survival Probabilities
#'
#' @param z A matrix or data frame of covariates. Each row corresponds to an
#'   observation, and each column corresponds to a predictor variable.
#' @param delta A numeric vector for the event indicator, where typically 1
#'   indicates an event occurred and 0 indicates censoring.
#' @param time A numeric vector of the observation times (either event or
#'   censoring time).
#' @param beta A numeric vector of regression coefficients, one for each
#'   column in `z`.
#' @param stratum An optional vector specifying the stratum for each observation.
#'   If `NULL` or missing, a single-stratum model is assumed.
#'
#' @return A matrix where each row corresponds to an individual and each column
#'   corresponds to an ordered event time. The value `S[i, j]` is the
#'   estimated survival probability for subject `i` at the `j`-th time point.
#'
#' @export
cal_surv_prob <- function(z, delta, time, beta, stratum) {
  if (missing(stratum)) {
    stratum <- as.matrix(rep(1, nrow(z)))
    colnames(stratum) <- "stratum"
    time_order <- order(time)
  } else {
    stratum <- as.matrix(match(stratum, unique(stratum)))
    time_order <- order(stratum, time)
  }
  
  delta <- as.matrix(delta[time_order])
  z <- as.matrix(z)[time_order, , drop = FALSE]
  stratum <- as.matrix(stratum[time_order, , drop = FALSE])
  beta <- as.matrix(beta)
  
  n.each_stratum <- as.numeric(table(stratum))
  
  diff <- ddloglik_S0(z, delta, beta, n.each_stratum)
  S0 <- diff$S0
  
  # Calculate hazard increments, handling potential division by zero
  hazard_increment <- rep(0, nrow(delta))
  if (any(S0 > 0)) {
    hazard_increment[S0 > 0] <- as.vector(delta)[S0 > 0] / S0[S0 > 0]
  }
  
  # Calculate the baseline cumulative hazard Lambda0 within each stratum
  Lambda0 <- ave(hazard_increment, as.vector(stratum), FUN = cumsum)
  
  tmax <- length(Lambda0)
  n <- nrow(z)
  S <- matrix(0, nrow = n, ncol = tmax)
  
  # Calculate individual survival probabilities S(t|Z) = exp(-Lambda0_k(t) * exp(Z'beta))
  risk_scores <- exp(z %*% beta)
  for (i in 1:tmax) {
    S[, i] <- exp(-Lambda0[i] * risk_scores)
  }
  return(S)
}


#' Calculate the Log-Likelihood for a Stratified Survival Model
#'
#' @param z A matrix or data frame of covariates. Each row corresponds to an
#'   observation, and each column corresponds to a predictor variable.
#' @param delta A numeric vector for the event indicator, where typically 1
#'   indicates an event occurred and 0 indicates censoring.
#' @param time A numeric vector of the observation times (either event or
#'   censoring time).
#' @param stratum An optional vector specifying the stratum for each observation.
#'   Can be a factor, character, or numeric vector. If `NULL` or missing, a
#'   single-stratum model is assumed.
#' @param beta A numeric vector of regression coefficients, one for each
#'   column in `z`.
#'
#' @return A single numeric value representing the calculated log-likelihood of
#'   the model.
#'
#' @export
loss_fn <- function(z, delta, time, stratum, beta) {
  if (missing(stratum)) {
    stratum <- as.matrix(rep(1, nrow(z)))
    colnames(stratum) <- "stratum"
    time_order <- order(time)
  } else {
    stratum <- as.matrix(match(stratum, unique(stratum)))
    time_order <- order(stratum, time)
  }
  
  delta <- as.matrix(delta[time_order])
  z <- as.matrix(z)[time_order, , drop = FALSE]
  stratum <- as.matrix(stratum[time_order, , drop = FALSE])
  beta <- as.matrix(beta)
  
  n.each_stratum <- as.numeric(table(stratum))
  
  diff <- loss_fn_cpp(z, delta, beta, n.each_stratum)
  loglik <- diff$loglik
  
  return(loglik)
}