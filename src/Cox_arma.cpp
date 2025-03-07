#include <RcppArmadillo.h>
#include "shared_functions.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// // [[Rcpp::export]]
// arma::vec rev_cumsum(const arma::vec& X){
//   int n = X.n_rows;
//   arma::vec rev_X = arma::flipud(X);
//   arma::vec cumsum_rev_X(n, arma::fill::zeros);
// 
//   double sum_to_i = 0;
//   for (int i = 0; i < n; ++i) {
//     sum_to_i += rev_X(i);
//     cumsum_rev_X(i) = sum_to_i;
//   }
// 
//   return arma::flipud(cumsum_rev_X);
// }

// [[Rcpp::export]]
arma::vec rev_cumsum(const arma::vec& X) {
  return arma::flipud(arma::cumsum(arma::flipud(X))); 
}


// [[Rcpp::export]]
List ddloglik(const arma::mat& Z, const arma::vec& delta, const arma::vec& beta){
  int p = beta.n_rows;
  int n = delta.n_rows;
  arma::vec loglik(1, arma::fill::zeros);
  arma::vec L1(p, arma::fill::zeros);
  arma::mat L2(p, p, arma::fill::zeros);
  arma::vec theta = Z * beta;
  arma::vec exp_theta = arma::exp(theta);
  arma::vec S0 = rev_cumsum(exp_theta);
  arma::mat S1(n, p, arma::fill::zeros);
  arma::mat S1_pre(n, p, arma::fill::zeros);

  for (int i = 0; i < p; ++i) {
    S1_pre.col(i) = Z.col(i) % exp_theta;
    S1.col(i) = rev_cumsum(S1_pre.col(i));
  }

  L1 = delta % (Z.each_col() - S1.each_col() / S0) * arma::ones<arma::mat>(p, 1);
  loglik(0) = arma::accu(delta % (theta - arma::log(S0)));

  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      arma::vec S2_i_colj = rev_cumsum(Z.col(j) % exp_theta % Z.col(i));
      arma::vec V_i_colj = (S2_i_colj / S0) - (S1.col(j) % S1.col(i) / arma::square(S0));
      L2(i, j) = arma::accu(delta % V_i_colj);
    }
  }

  return List::create(
    Named("loglik") = loglik,
    Named("L1") = L1,
    Named("L2") = L2,
    Named("S0") = S0
  );
}


// [[Rcpp::export]]
List ddloglik_S0(const arma::mat& Z, const arma::vec& delta, const arma::vec& beta){
  arma::vec theta = Z * beta;
  arma::vec exp_theta = arma::exp(theta);
  arma::vec S0 = rev_cumsum(exp_theta);
  
  return List::create(
    Named("S0") = S0
  );
}


