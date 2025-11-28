#' Calculate Survival Probabilities
#'
#' @description
#' Computes individual survival probabilities from a fitted linear predictor
#' \code{z\%*\%beta} using a stratified Breslow-type baseline hazard estimate.
#'
#' @param z A numeric matrix (or data frame coercible to matrix) of covariates.
#'   Each row is an observation and each column a predictor.
#' @param delta A numeric vector of event indicators (1 = event, 0 = censored).
#' @param time A numeric vector of observed times (event or censoring).
#' @param beta A numeric vector of regression coefficients with length equal to
#'   the number of columns in \code{z}.
#' @param stratum An optional vector specifying the stratum for each observation. 
#'   If missing, a single-stratum model is assumed.
#'
#' @details
#' Inputs are internally sorted by \code{stratum} and \code{time}. Within each
#' stratum, a baseline hazard increment is computed as \code{delta/S0}, where
#' \code{S0} is the risk set sum returned by \code{ddloglik_S0}. The stratified
#' baseline cumulative hazard \code{Lambda0} is then formed by a cumulative sum
#' within stratum, and individual survival curves are computed as
#' \code{S(t) = exp(-Lambda0(t) * exp(z \%*\% beta))}.
#'
#' @return
#' A numeric matrix of survival probabilities with \code{nrow(z)} rows and
#' \code{length(time)} columns. Rows correspond to observations; columns are in
#' the internal sorted order of \code{(stratum, time)} (i.e., not collapsed to
#' unique event times). Entry \code{S[i, j]} is the estimated survival
#' probability for subject \code{i} evaluated at the \code{j}-th sorted time
#' point.
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


#' Calculate the Log-Partial Likelihood for a Stratified Cox Model
#'
#' @description
#' Computes the stratified Cox partial log-likelihood for given covariates,
#' event indicators, times, and coefficients.
#'
#' @param z A numeric matrix (or data frame coercible to matrix) of covariates.
#'   Each row is an observation and each column a predictor.
#' @param delta A numeric vector of event indicators (1 = event, 0 = censored).
#' @param time A numeric vector of observed times (event or censoring).
#' @param stratum An optional vector specifying the stratum for each observation
#'   (factor/character/numeric). If missing, a single-stratum model is assumed.
#' @param beta A numeric vector of regression coefficients with length equal to
#'   the number of columns in \code{z}.
#'
#' @details
#' Inputs are internally sorted by \code{stratum} and \code{time}. The function
#' evaluates the stratified Cox partial log-likelihood using the supplied \code{z}, 
#' \code{delta}, \code{beta}, and the stratum sizes.
#'
#' @return
#' A single numeric value giving the stratified Cox partial log-likelihood.
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