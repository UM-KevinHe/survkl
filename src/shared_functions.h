// shared_functions.h
#ifndef SHARED_FUNCTIONS_H
#define SHARED_FUNCTIONS_H

#include <RcppArmadillo.h>

// Declare the function
arma::vec rev_cumsum(const arma::vec& X);

arma::vec calculateDeltaTilde(const arma::vec& event, const arma::vec& time, const arma::vec& theta_tilde);

Rcpp::List calculateRiskAndUpdateLoss(arma::vec& haz, const arma::vec& d, double Loss, arma::vec r);

Rcpp::List ddloglik_KL_RS_score(const arma::mat& Z, const arma::vec& delta, arma::vec& beta, const arma::vec& theta_tilde, const double &eta);

#endif // SHARED_FUNCTIONS_H
