#include <RcppArmadillo.h>
#include "shared_functions.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
double ddloglik_KL_RS_test(double z){
  return z;
}

// [[Rcpp::export]]
List loss_fn_cpp(const arma::mat& Z, const arma::vec& delta, arma::vec& beta) {
  int n = delta.n_rows;
  double loglik = 0;
  
  arma::vec theta = Z * beta;
  arma::vec exp_theta = arma::exp(theta);
  arma::vec S0 = rev_cumsum(exp_theta);

  loglik = arma::accu(delta % (theta - arma::log(S0)));
  
  //print L2 dimensions:
  // Rcpp::Rcout << "L2 dimensions: " << L2.n_rows << " x " << L2.n_cols << std::endl;
  
  return List::create(Named("loglik") = loglik);
}


// [[Rcpp::export]]
List ddloglik_KL_RS_score(const arma::mat& Z, const arma::vec& delta, arma::vec& beta, const arma::vec& theta_tilde, const double &eta) {
  int p = beta.n_rows;
  int n = delta.n_rows;
  double loglik = 0;
  arma::vec L1 = arma::zeros<arma::vec>(p);
  arma::vec theta = Z * beta;
  arma::vec exp_theta = arma::exp(theta);
  arma::vec S0 = rev_cumsum(exp_theta);
  arma::mat S1 = arma::zeros<arma::mat>(n, p);

  arma::vec exp_theta_tilde = arma::exp(theta_tilde);
  arma::vec S0_tilde = rev_cumsum(exp_theta_tilde);
  arma::mat S1_tilde = arma::zeros<arma::mat>(n, p);
  
  for (int i = 0; i < p; i++) {
    arma::vec S1_pre = Z.col(i) % exp_theta;
    S1.col(i) = rev_cumsum(S1_pre);
    
    arma::vec S1_pre_tilde = Z.col(i) % exp_theta_tilde;
    S1_tilde.col(i) = rev_cumsum(S1_pre_tilde);
  }
  
  // Perform the operation in steps for clarity
  arma::mat temp_1(n, p, arma::fill::zeros); 
  
  for (int i = 0; i < p; i++) {
    temp_1.col(i) = delta % (Z.col(i) + eta * (S1_tilde.col(i) / S0_tilde) - (1 + eta) * (S1.col(i) / S0));
  }
  
  L1 = arma::sum(temp_1, 0).t();
  loglik = arma::accu(delta % (theta - arma::log(S0)));
  
  
  return List::create(Named("loglik") = loglik,
                      Named("L1") = L1);
}



static inline void score_hess_kl(
  const arma::mat& Z,
  const arma::vec& delta,
  const arma::vec& beta,
  const arma::vec& exp_theta_tilde,   // precomputed exp(theta_tilde)
  const arma::vec& S0_tilde,          // precomputed rev_cumsum(exp_theta_tilde)
  const arma::mat& S1_tilde,          // precomputed rev_cumsum(Z.col(k) % exp_theta_tilde), n x p
  const double eta,
  arma::vec& L1,
  arma::mat& L2
){
  const int n = Z.n_rows;
  const int p = Z.n_cols;

  arma::vec theta = Z * beta;
  arma::vec exp_theta = arma::exp(theta);
  arma::vec S0 = rev_cumsum(exp_theta);

  // Build S1 (n x p): rev_cumsum(Z.col(j) % exp_theta) colwise
  arma::mat S1(n, p, arma::fill::zeros);
  for (int j = 0; j < p; ++j) {
    S1.col(j) = rev_cumsum(Z.col(j) % exp_theta);
  }

  // Score
  // temp_1.col(j) = delta % ( Z.col(j) + eta * (S1_tilde.col(j)/S0_tilde) - (1+eta) * (S1.col(j)/S0) )
  arma::mat temp_1(n, p, arma::fill::zeros);
  for (int j = 0; j < p; ++j) {
    temp_1.col(j) = delta % ( Z.col(j) 
                              + eta * (S1_tilde.col(j) / S0_tilde)
                              - (1.0 + eta) * (S1.col(j) / S0) );
  }
  L1 = arma::sum(temp_1, 0).t();

  L2.zeros();
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      arma::vec S2_jk = rev_cumsum( (Z.col(j) % Z.col(k)) % exp_theta );
      arma::vec V_jk  = (S2_jk / S0) - (S1.col(j) % S1.col(k)) / arma::square(S0);
      L2(j,k) = (1.0 + eta) * arma::accu( delta % V_jk );
    }
  }

  // Small ridge for numerical stability (doesn't change the target)
  L2.diag() += 1e-8;
}

// [[Rcpp::export]]
arma::vec KL_Cox_Estimate_cpp(
  const arma::mat& Z, 
  const arma::vec& delta, 
  const arma::vec& theta_tilde,  // external linear predictor (risk score), length n
  const double eta,
  const double tol = 1.0e-7,
  const int maxit = 50
){
  const int p = Z.n_cols;
  arma::vec beta = arma::zeros<arma::vec>(p);

  // Precompute external part once
  arma::vec exp_theta_tilde = arma::exp(theta_tilde);
  arma::vec S0_tilde        = rev_cumsum(exp_theta_tilde);
  arma::mat S1_tilde(Z.n_rows, p, arma::fill::zeros);
  for (int j = 0; j < p; ++j) {
    S1_tilde.col(j) = rev_cumsum(Z.col(j) % exp_theta_tilde);
  }

  // Newton
  for (int it = 0; it < maxit; ++it) {
    arma::vec L1(p, arma::fill::zeros);
    arma::mat L2(p, p, arma::fill::zeros);
    score_hess_kl(Z, delta, beta, exp_theta_tilde, S0_tilde, S1_tilde, eta, L1, L2);

    // Solve L2 * step = L1
    arma::vec step;
    bool ok = arma::solve(step, L2, L1, arma::solve_opts::likely_sympd);
    if (!ok) {
      arma::mat L2r = L2;
      L2r.diag() += 1e-4;
      ok = arma::solve(step, L2r, L1, arma::solve_opts::likely_sympd);
      if (!ok) Rcpp::stop("Hessian solve failed");
    }

    beta += step;
    if (arma::max(arma::abs(step)) < tol) break;
  }

  return beta;
}

// [[Rcpp::export]]
arma::vec KL_Cox_Estimate_cpp_ridge(
  const arma::mat& Z, 
  const arma::vec& delta, 
  const arma::vec& theta_tilde,
  const double eta,
  const double lambda,           // ridge strength (>=0)
  const double tol = 1.0e-7,
  const int maxit = 50
){
  const int p = Z.n_cols;
  arma::vec beta = arma::zeros<arma::vec>(p);

  arma::vec exp_theta_tilde = arma::exp(theta_tilde);
  arma::vec S0_tilde        = rev_cumsum(exp_theta_tilde);
  arma::mat S1_tilde(Z.n_rows, p, arma::fill::zeros);
  for (int j = 0; j < p; ++j) {
    S1_tilde.col(j) = rev_cumsum(Z.col(j) % exp_theta_tilde);
  }

  for (int it = 0; it < maxit; ++it) {
    arma::vec L1(p, arma::fill::zeros);
    arma::mat L2(p, p, arma::fill::zeros);
    score_hess_kl(Z, delta, beta, exp_theta_tilde, S0_tilde, S1_tilde, eta, L1, L2);

    arma::mat H = L2;
    H.diag() += lambda;
    arma::vec rhs = L1 - lambda * beta;

    arma::vec step;
    bool ok = arma::solve(step, H, rhs, arma::solve_opts::likely_sympd);
    if (!ok) {
      H.diag() += 1e-4;
      ok = arma::solve(step, H, rhs, arma::solve_opts::likely_sympd);
      if (!ok) Rcpp::stop("Hessian solve failed (even with ridge).");
    }

    beta += step;
    if (arma::max(arma::abs(step)) < tol) break;
  }
  return beta;
}