#' Simulate Low-Dimensional Survival Data for Internal and External Cohorts
#'
#' @description
#' Generates simulated low-dimensional survival datasets for internal (training/testing) 
#' and external cohorts with varying heterogeneity levels.
#'
#' @param n_int Number of subjects in the internal training set (default = 200).
#' @param n_test Number of subjects in the internal test set (default = 1000).
#' @param n_ext Number of subjects in the external dataset (default = 1000).
#' @param beta_true True regression coefficients.
#' @param int_cens_target Target censoring rate for internal data.
#' @param ext_cens_target Target censoring rate for external data.
#' @param lambda0,nu0 Parameters of the Weibull baseline hazard function 
#'   \eqn{h_0(t) = \lambda_0 \nu_0 t^{\nu_0 - 1}}.
#' @param qualities Character; external data quality level, one of 
#'   \code{"Good"}, \code{"Fair"}, or \code{"Poor"}.
#' @param pL_map_external latent mixture for external, indicate the heterogeneity between external and internal cohorts.
#' @param seed Random seed for reproducibility.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{external}{External dataset.}
#'   \item{internal_train}{Internal training dataset.}
#'   \item{internal_test}{Internal test dataset.}
#' }
#'
#' @examples
#' data_list <- sim_lowdim()
#' str(data_list)
#'
sim_lowdim <- function(n_int = 200, n_test = 1000, n_ext = 1000, 
                       beta_true = c(0.3, -0.3,  0.3, -0.3,  0.3, -0.3),
                       int_cens_target = 0.3, ext_cens_target = 0.5, 
                       lambda0 =  1, nu0 = 2,
                       qualities = c("Good","Fair","Poor"),
                       pL_map_external = c(Good = 1.0, Fair = 0.5, Poor = 0.0), 
                       seed = 1){
  
  set.seed(seed)
  qualities <- match.arg(qualities, c("Good","Fair","Poor"), several.ok=FALSE)
  
  names(beta_true) <- paste0("Z", 1:6)
  pL_ext <- unname(pL_map_external[qualities])
  df_ext <- sim(n_ext, ext_cens_target, beta_true = beta_true, pL_group = pL_ext)
  
  pL_internal <- 1.0
  df_tr <- sim(n_int, int_cens_target, beta_true = beta_true, pL_group = pL_internal)
  df_te <- sim(n_test, int_cens_target, beta_true = beta_true, pL_group = pL_internal)

  return(list(external = df_ext, 
              internal_train=df_tr, 
              internal_test=df_te))
}


## helper functions
simulate_survival_weibull <- function(eta, lambda = lambda0, nu = nu0) {
  U <- runif(length(eta))
  (-log(U) / (lambda * exp(eta)))^(1/nu)
}


sim <- function(n, target_cens, beta_true, pL_group = 1.0) {
  Z1Z2 <- simu_z_AR1(n, 2, rho=0.5)
  colnames(Z1Z2) <- c("Z1","Z2")
  Z3   <- rbinom(n, 1, 0.5)
  Z4   <- rbinom(n, 1, 0.5)
  L    <- rbinom(n, 1, pL_group)  # a length-n vector of latent group indicator (0/1)
  Z5   <- rnorm(n, mean= 2*L, sd=1)  # Z5 & Z6 from a mixture of two populations
  Z6   <- rnorm(n, mean=-2*L, sd=1)
  X    <- cbind(Z1Z2, Z3=Z3, Z4=Z4, Z5=Z5, Z6=Z6)
  eta  <- drop(as.matrix(X) %*% beta_true)
  Tt   <- simulate_survival_weibull(eta)
  ub   <- tune_censoring_uniform(Tt, target=target_cens)
  Cc   <- runif(n, 0, ub)
  status <- as.integer(Tt <= Cc)
  time   <- ifelse(status==1, Tt, Cc)
  order_cols(data.frame(X, status=status, time=time), z_names=paste0("Z",1:6))
}


AR1_corr <- function(rho, m) {
  if (m == 1) return(matrix(1,1,1))
  out <- outer(1:m, 1:m, function(i,j) rho^abs(i-j))
  diag(out) <- 1
  out
}

simu_z_AR1 <- function(n, p_cont, rho=0.5) {
  mvtnorm::rmvnorm(n, mean=rep(0,p_cont), sigma=AR1_corr(rho, p_cont))
}


tune_censoring_uniform <- function(event_time, target=0.5, max_iter=30) {
  q <- stats::quantile(event_time, probs=c(0.2, 0.5, 0.95), type=7, na.rm=TRUE, names=FALSE)
  ub_min <- if (is.finite(q[1]) && q[1] > 0) q[1] else 1e-3
  ub_max <- if (is.finite(q[3]) && q[3] > 0) q[3] else max(1, mean(event_time, na.rm=TRUE))
  ub     <- if (is.finite(q[2]) && q[2] > 0) q[2] else mean(event_time[event_time>0], na.rm=TRUE)
  if (!is.finite(ub) || ub <= 0) ub <- 1
  ub <- min(max(ub, ub_min), ub_max)
  best_ub <- ub; best_diff <- Inf
  for (i in 1:max_iter) {
    C <- runif(length(event_time), 0, ub)
    delta <- as.integer(event_time <= C)
    rate  <- mean(1 - delta) # censoring rate
    d <- abs(rate - target)
    if (d < best_diff) { best_diff <- d; best_ub <- ub }
    ub <- if (rate < target) ub*0.85 else ub*1.15
    ub <- min(max(ub, ub_min), ub_max)
  }
  best_ub
}


order_cols <- function(df, z_names) {
  keep <- intersect(z_names, names(df))
  miss <- setdiff(z_names, keep)
  if (length(miss)) for (m in miss) df[[m]] <- 0
  df <- df[, c(z_names, "status", "time"), drop=FALSE]
  rownames(df) <- NULL
  df
}