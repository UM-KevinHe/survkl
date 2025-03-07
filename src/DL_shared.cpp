#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "shared_functions.h"


// [[Rcpp::export]]
arma::vec calculateDeltaTilde(const arma::vec& event, const arma::vec& time, const arma::vec& theta_tilde) {
  int n = event.size();
  arma::vec failure_time = time.elem(arma::find(event == 1));
  int m = failure_time.size();
  
  arma::vec stabilized_exp_theta_tilde = exp(theta_tilde);
  
  arma::vec cumsum_reversed_stabilized = arma::flipud(arma::cumsum(arma::flipud(stabilized_exp_theta_tilde)));
  arma::vec denominator = 1 / cumsum_reversed_stabilized;
  arma::vec denominator_kfailures = denominator.elem(arma::find(event == 1));
  
  arma::vec delta_tilde = arma::zeros<arma::vec>(n);
  
  for (int k = 0; k < m; ++k) {
    arma::vec Y_i = arma::conv_to<arma::vec>::from(time >= failure_time[k]);
    arma::vec numerator = Y_i % stabilized_exp_theta_tilde;
    delta_tilde += numerator * denominator_kfailures[k];
  }
  
  return delta_tilde;
}

// [[Rcpp::export]]
List calculateRiskAndUpdateLoss(arma::vec& eta, const arma::vec& d, 
                                double Loss, arma::vec r) {
  int n = eta.n_elem;
  arma::vec haz = arma::zeros<arma::vec>(n), rsk = arma::zeros<arma::vec>(n), h = arma::zeros<arma::vec>(n);
  double v = 1.0, s = 0.0;
  
  // Calculate haz, risk
  Rcpp::Rcout << "inside haz: " << accu(haz) << std::endl;
  rsk(n-1) = haz(n-1);
  for (int i = n-2; i >= 0; i--) {
    rsk(i) = rsk(i+1) + haz(i);
  }
  for (int i = 0; i < n; i++) {
    Loss += d(i) * eta(i) - d(i) * log(rsk(i));
  }
  // Approximate L:
  h(0) = d(0) / rsk(0);
  for (int i = 1; i < n; i++) {
    h(i) = h(i-1) + d(i) / rsk(i);
  }
  for (int i = 0; i < n; i++) {
    h(i) = h(i) * haz(i);
    s = d(i) - h(i);
    if (h(i) == 0) r(i) = 0;
    else r(i) = s / v;
  }
  
  Rcpp::Rcout << "inside Loss: " << Loss << std::endl;
  Rcpp::Rcout << "inside r: " << r(0) << std::endl;
  
  return List::create(Named("Loss") = Loss, Named("r") = r);
}


// [[Rcpp::export]]
List calculateRiskAndUpdateLoss2(arma::vec& eta, arma::vec& d, arma::vec& deltaTilde, double& Loss, double& gamma) {

  int n = eta.size();
  arma::vec haz(n), rsk(n), h(n), r(n);
  double v = 1.0, s = 0.0;
  
  // Calculate haz, risk
  haz = exp(eta);
  rsk(n-1) = haz(n-1);
  for (int i = n-2; i >= 0; i--) {
    rsk(i) = rsk(i+1) + haz(i);
  }
  for (int i = 0; i < n; i++) {
    Loss += d(i) * eta(i) - d(i) * log(rsk(i));
  }
  // Approximate L:
  h(0) = d(0) / rsk(0);
  for (int i = 1; i < n; i++) {
    h(i) = h(i-1) + d(i) / rsk(i);
  }
  for (int i = 0; i < n; i++) {
    h(i) = h(i) * haz(i);
    s = (d(i)+gamma*deltaTilde(i))/(1+gamma) - h(i);
    if (h(i) == 0) r(i) = 0;
    else r(i) = s / v;
  }
  
  return List::create(Named("Loss") = Loss, Named("r") = r);
  
}

// // [[Rcpp::export]]
// List DL_KL_score(const arma::mat& Z, const arma::vec& delta, arma::vec& beta, const arma::vec& theta_tilde, const arma::vec& delta_tilde, double &eta) {
//   int p = beta.n_rows;
//   int n = delta.n_rows;
//   double loglik = 0;
//   arma::vec L1 = arma::zeros<arma::vec>(p);
//   arma::mat L2 = arma::zeros<arma::mat>(p, p);
//   arma::vec theta = Z * beta;
//   arma::vec exp_theta = arma::exp(theta);
//   arma::vec S0 = rev_cumsum(exp_theta);
//   arma::mat S1 = arma::zeros<arma::mat>(n, p);
// 
//   for (int i = 0; i < p; i++) {
//     arma::vec S1_pre = Z.col(i) % exp_theta;
//     S1.col(i) = rev_cumsum(S1_pre);
//   }
// 
//   // Perform the operation in steps for clarity
//   arma::mat temp_1(n, p, arma::fill::zeros); 
// 
//   for (int i = 0; i < p; i++) {
//     temp_1.col(i) = delta % (Z.col(i) + eta * (S1_tilde.col(i) / S0_tilde) - (1 + eta) * (S1.col(i) / S0));
//   }
// 
//   L1 = arma::sum(temp_1, 0).t();
//   loglik = arma::accu(delta % (theta - arma::log(S0)));
//   ///
// 
//   return List::create(Named("loglik") = loglik,
//                       Named("L1") = L1);
// }