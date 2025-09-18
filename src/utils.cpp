#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "utils.h"

using namespace Rcpp;


/*
 * Computes the reverse cumulative sum of a numeric vector.
 * 
 * Parameters:
 *   X - Input numeric vector (length n)
 * 
 * Example:
 *   X = (1, 2, 3)
 *   rev_cumsum(X) = (6, 5, 3)
 */
// [[Rcpp::export]]
arma::vec rev_cumsum(const arma::vec& X) {
  return arma::flipud(arma::cumsum(arma::flipud(X))); 
}


/*
 * Computes the log-partial likelihood for the stratified Cox model.
 *
 * Parameters:
 *   lp              - Linear predictor (length n), i.e., Z * beta
 *   delta           - Event indicator vector (length n), 1 = event, 0 = censored
 *   
 *   n_each_stratum  - Vector indicating number of observations per stratum
 *
 * Assumes:
 *   - Data is sorted by stratum first, then by increasing time.
 *
 * Returns:
 *   The log-partial likelihood value (scalar)
 */
// [[Rcpp::export]]
double pl_cal_theta(const arma::vec& lp,
                    const arma::vec& delta,
                    const arma::vec& n_each_stratum) {

  const arma::uword S = n_each_stratum.n_elem;
  arma::uvec ind_start(S);
  ind_start(0) = 0;
  for (arma::uword s = 1; s < S; ++s) {
    ind_start(s) = ind_start(s - 1) + n_each_stratum(s - 1);
  }

  double loglik = 0.0;

  for (arma::uword s = 0; s < S; ++s) {
    arma::uword start = ind_start(s);
    arma::uword len = n_each_stratum(s);
    arma::uword end = start + len - 1;

    arma::vec lp_s = lp.subvec(start, end);
    arma::vec delta_s = delta.subvec(start, end);
    arma::vec exp_lp_s = arma::exp(lp_s);

    arma::vec S0_s = rev_cumsum(exp_lp_s);

    for (arma::uword i = 0; i < len; ++i) {
      if (delta_s(i) > 0) {
        loglik += lp_s(i) - std::log(S0_s(i));
      }
    }
  }
  return loglik;
}


/* * Computes the delta_tilde vector for the Cox-KL model.
 * This vector is used to calculate the predicted event indicator vector, defined using external risk scores.
 *
 * Parameters:
 *   event         - Event indicator vector (length n), 1 = event, 0 = censored
 *   time          - Time-to-event vector (length n)
 *   RS   - Risk score based on external covariates (length n), i.e., Z * beta_tilde
 *   n_each_stratum - Vector indicating the number of observations in each stratum
 *
 * Returns:
 *   A vector of adjusted risk scores (delta_tilde)
 */
// [[Rcpp::export]]
arma::vec calculateDeltaTilde(const arma::vec& event, 
                              const arma::vec& time, 
                              const arma::vec& RS,
                              const arma::vec& n_each_stratum) {
  const arma::uword n = event.n_elem;
  arma::vec delta_tilde(n, arma::fill::zeros);    
  arma::vec exp_theta_tilde = arma::exp(RS);
  
  const arma::uword S = n_each_stratum.n_elem;
  arma::uvec ind_start(S); 
  ind_start(0) = 0;
  for (arma::uword i = 1; i < S; ++i) {
    ind_start(i) = ind_start(i - 1) + n_each_stratum(i - 1);
  }
  
  for (arma::uword s = 0; s < S; ++s){
    arma::uword start = ind_start(s);
    arma::uword len   = n_each_stratum(s);
    if (len == 0) continue;
    
    const arma::uword end = start + len - 1;
    
    arma::vec event_s = event.subvec(start, end);
    arma::vec time_s = time.subvec(start, end);
    arma::vec exp_theta_tilde_s = exp_theta_tilde.subvec(start, end);
    
    arma::vec denom_s = rev_cumsum(exp_theta_tilde_s);
    
    arma::uvec failure_idx_s = arma::find(event_s == 1);
    if (failure_idx_s.is_empty()) continue;
    
    for (arma::uword k = 0; k < failure_idx_s.n_elem; ++k) {
      const arma::uword f = failure_idx_s(k);
      arma::uvec at_risk = arma::find(time_s >= time_s(f));
      if (at_risk.is_empty()) continue;
      
      const double denom = denom_s(f);
      if (denom <= 0.0) continue;
      
      arma::uvec at_risk_global = at_risk + start;
      delta_tilde.elem(at_risk_global) += exp_theta_tilde_s.elem(at_risk) / denom;
    }
  }
  
  
  return delta_tilde;
}


/*
 * Computes the Cox partial log-likelihood given covariates and parameters.
 *
 * Parameters:
 *   Z              - Covariate matrix (n × p)
 *   delta          - Event indicator vector (length n, 1 = event, 0 = censored)
 *   beta           - Regression coefficient vector (p × 1)
 *   n_each_stratum - Number of observations in each stratum (length = # strata)
 *
 * Returns:
 *   A list with:
 *     - loglik : value of the Cox partial log-likelihood
 *
 * Notes:
 *   - The calculation is stratified
 */
// [[Rcpp::export]]
List loss_fn_cpp(const arma::mat& Z, 
                 const arma::vec& delta,
                 arma::vec& beta,
                 const arma::vec& n_each_stratum) {
  
  const arma::uword S = n_each_stratum.n_elem;
  arma::vec theta = Z * beta;
  arma::vec exp_theta = arma::exp(theta);
  double loglik = 0.0;
  
  // Compute start index for each stratum
  arma::vec ind_start(S);
  ind_start(0) = 0;
  for (arma::uword i = 1; i < S; ++i) {
    ind_start(i) = ind_start(i - 1) + n_each_stratum(i - 1);
  }
  
  // Loop over strata
  for (arma::uword j = 0; j < S; ++j) {
    arma::uword start = ind_start(j);
    arma::uword len = n_each_stratum(j);
    arma::uword end = start + len - 1;

    arma::vec theta_s = theta.subvec(start, end);
    arma::vec exp_theta_s = exp_theta.subvec(start, end);
    arma::vec delta_s = delta.subvec(start, end);
    
    arma::vec S0_s = rev_cumsum(exp_theta_s);
    loglik += arma::accu(delta_s % (theta_s - arma::log(S0_s)));
  }
  
  return List::create(Named("loglik") = loglik);
}



/*
 * Computes S0 for the stratified Cox models.
 *
 * Definition:
 *   For each observation i within a stratum s, 
 *   S0_i = sum_{l in stratum s : t_l >= t_i} exp(Z * beta).
 *
 * Parameters:
 *   Z               - Covariate matrix (n × p)
 *   delta           - Event indicator vector (length n)  [not used here; included for interface symmetry]
 *   beta            - Regression coefficients (p × 1)
 *   n_each_stratum  - Number of observations in each stratum (length = # strata)
 *
 * Returns:
 *   A list with:
 *     - S0 : vector of length n containing the stratum-wise reverse cumulative sums of exp(theta).
 */
// [[Rcpp::export]]
List ddloglik_S0(const arma::mat& Z, 
                 const arma::vec& delta, 
                 const arma::vec& beta, 
                 const arma::vec& n_each_stratum) {
  const arma::uword S = n_each_stratum.n_elem;
  arma::vec theta = Z * beta;
  arma::vec exp_theta = arma::exp(theta);
  arma::vec S0(Z.n_rows);  // Output
  
  // Compute start index for each stratum
  arma::vec ind_start(S);
  ind_start(0) = 0;
  for (arma::uword i = 1; i < S; ++i) {
    ind_start(i) = ind_start(i - 1) + n_each_stratum(i - 1);
  }

  for (arma::uword j = 0; j < S; ++j) {
    arma::uword start = ind_start(j);
    arma::uword len = n_each_stratum(j);
    arma::uword end = start + len - 1;
    
    arma::vec exp_theta_s = exp_theta.subvec(start, end);
    arma::vec S0_s = rev_cumsum(exp_theta_s);
    S0.subvec(start, end) = S0_s;
  }
  
  return List::create(Named("S0") = S0);
}



/*
 * Computes the concordance index (C-index) for a single stratum of a Cox model.
 *
 * Parameters:
 *   time   - Observed time-to-event vector (length n), sorted in increasing order
 *   xbeta  - Linear predictor vector (length n), i.e., Z * beta
 *   delta  - Event indicator vector (length n), 1 = event, 0 = censored
 *
 * Assumes:
 *   - Input vectors are sorted by increasing time.
 *   - Ties in event time or predicted risk are handled with 0.5 weighting.
 *   - Only pairs where the earlier individual experienced an event contribute.
 *
 * Returns:
 *   A named list with:
 *     - numer: number of concordant pairs (possibly fractional)
 *     - denom: number of comparable pairs
 */

 // [[Rcpp::export]]
Rcpp::List cox_c_index(const arma::vec& time,
                       const arma::vec& xbeta,
                       const arma::vec& delta) {
  int N = time.n_elem;
  double factor = (N > 100000) ? std::pow(10.0, -8) : 1.0;

  double N_pair = 0.0;  // Total number of comparable pairs
  double N_con  = 0.0;  // Total number of concordant pairs

  for (int i = 0; i < N; ++i) {
    if (delta(i) == 1) {
      N_pair += (N - i - 1) * factor;

      for (int j = i + 1; j < N; ++j) {
        double obs_diff = std::abs(time(i) - time(j));
        double risk_diff = std::abs(xbeta(i) - xbeta(j));

        if (obs_diff == 0) {
          if (delta(j) == 1) {
            N_con += 0.5 * factor;
          } else {
            if (risk_diff < 1e-4) {
              N_con += 0.5 * factor;
            } else if (xbeta(i) > xbeta(j)) {
              N_con += factor;
            }
          }
        } else {
          if (risk_diff < 1e-4) {
            N_con += 0.5 * factor;
          } else if (xbeta(i) > xbeta(j)) {
            N_con += factor;
          }
        }
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("numer") = N_con,
    Rcpp::Named("denom") = N_pair
  );
}



// [[Rcpp::export]]
double maxgrad(arma::mat &x, arma::vec &r, arma::vec &K, arma::vec &m){ // "K": a vector contains the start index of each group; "m": group.multiplier
  int J = K.n_elem - 1;  //number of penalized group
  double z_max = 0, z;
  for (int g = 0; g < J; g++){
    int Kg = K(g + 1) - K(g); //number of features in group g
    arma::vec Z(Kg);
    for (int j = K(g); j < K(g + 1); j++) {
      arma::vec x_tmp = x.col(j);
      Z(j - K(g)) = dot(x_tmp, r);
    }
    z = arma::norm(Z)/m(g);
    if (z > z_max){
      z_max = z;
    }
  }
  return(z_max);
}