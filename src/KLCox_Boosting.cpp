#include <RcppArmadillo.h>
#include <math.h>
#include "shared_functions.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace Rcpp;

//' An Rcpp function that calculates the first and second order derivate and partial likelihood function.
//' @param N   The sample size
//' @param delta   Event indicator
//' @param z   Covariate matrix
//' @param beta   The coefficient
//' @param offset   The whole spanned matrix
//' 
//' @return the first and second derivate and partial likelihood
// [[Rcpp::export]]
List ddloglik_md(int n, arma::colvec &delta, arma::mat &z, arma::colvec &beta, arma::mat &offset){
  int col   = beta.n_rows;     // row number of beta

  arma::mat s0    = arma::flipud(arma::cumsum(arma::flipud(exp(z*beta+offset))));
  arma::mat s00   = arma::repmat(s0,1,col);      // span the matrix in block structure
  arma::mat s0_2  = arma::repmat(s0%s0,1,col);   // for calculating the hessian matrix
  arma::mat delta0= arma::repmat(delta,1,col);   // span delta to calculating the first derivative
  arma::mat bb    = arma::repmat(exp(z*beta+offset),1,col);  
  arma::mat s     = z%bb;                   
  arma::mat s1    = arma::flipud(arma::cumsum(arma::flipud(s)));
  arma::mat l1    = (z-s1/s00)%delta0;        
  arma::rowvec L1 = arma::sum(l1,0);            // get the first derivative
  arma::mat lambda= arma::cumsum(delta/s0);     // get the lambda
  arma::mat partial_likelihood = delta%((z*beta+offset)-log(s0)); // get the partial likelihood (according to definination of partial likelihood of cox model)
  // arma::mat s2;         
  // arma::mat l2;         
  // arma::mat L2= arma::zeros<arma::mat>(col,col);// initilize Hessian matrix
  // // Using a loop to calculate the Hessian Matrix
  // for(int i = 0; i < col; i++)
  // { s2      = s%(arma::repmat(z.col(i),1,col));
  //   s2      = arma::flipud(arma::cumsum(arma::flipud(s2)));
  //   arma::mat s11 = arma::repmat(s1.col(i),1,col);
  //   l2      = delta0%(s2/s00-s1%s11/(s0_2));
  //   L2.col(i) = sum(l2,0).t();
  // }
  // return the result we need for calculation.
  List result;
  result["L1"]  = L1;
  result["s0"]  = s0;
  result["s00"]  = s00;
  result["s0_2"] = s0_2;
  result["delta0"] =delta0;
  result["bb"]    = bb;
  // result["L2"]  = L2;
  result["lambda"] = lambda;
  result["partial_likelihood"] = partial_likelihood;
  return result;
}


// [[Rcpp::export]]
List ddloglik_md2(arma::vec &delta, arma::mat &z, arma::vec &beta){
  int n   = delta.n_rows;
  int p   = beta.n_rows;     // row number of beta
  
  arma::mat L         = arma::zeros<arma::mat>(n,p);
  arma::rowvec L1     = arma::zeros<arma::rowvec>(p);

  arma::vec S0     = arma::zeros<arma::vec>(n);
  arma::mat S1;
  arma::mat partial_likelihood;
  //test
  arma::mat zbeta;
  arma::vec S00;    //for each strata
  arma::mat    S10;    //for each strata
  arma::vec Lambda(n);
  arma::mat    matLambda;
  arma::mat    matdelta;

  //main function
  zbeta  = exp(z * beta);                                 //exp(beta * z)
  S00    = arma::reverse(arma::cumsum(arma::reverse(zbeta)));  
  S10    = z % repmat(zbeta,1,p);
  Lambda = cumsum(delta/S00);
  matLambda = repmat(Lambda,1,p);
  matdelta  = repmat(delta,1,p);
  L     = matdelta%z - S10%matLambda;
  L1 = sum(L,0);
  arma::colvec L11= L1.t();                               // get L1 norm

  partial_likelihood = delta%((z*beta)-log(S00));
  double likelihood = arma::accu(partial_likelihood);

  List result;
  result["L1"]  = L11;
  //result["L"]  = L;
  //result["S00"]  = S00;
  //result["S10"]  = S10;
  //result["Lambda"]  = Lambda;
  //result["matLambda"]  = matLambda;
  //result["matdelta"]  = matdelta;
  result["likelihood"] = likelihood;
  return result;
}




// [[Rcpp::export]]
List klcox_boosting(const arma::mat& z, const arma::vec& delta,
                    const arma::vec& theta_tilde, const double &eta, 
                    const double& rate, const double& tol, const int& maxit) {
  int n   = z.n_rows;
  int p   = z.n_cols;     // row number of beta
  
  arma::vec beta = arma::zeros<arma::vec>(p);
  int key = 0;
  bool converge = false;
  double loglik = 0;
  arma::vec GD;
  std::vector<double> likelihood_all;
  std::vector<int> j_star_all;
  
  while (!converge) {
    key++;
    List result = ddloglik_KL_RS_score(z, delta, beta, theta_tilde, eta);
    arma::vec update_all = result["L1"];
    loglik = result["loglik"];

    GD = square(update_all);
    int j_star = index_max(GD);
    j_star_all.push_back(j_star + 1); // Adjust for 1-based indexing in R

    beta(j_star) += rate * (update_all(j_star) > 0 ? 1 : -1);
    likelihood_all.push_back(loglik);

    if (key >= 10) {
      double llk_diff = std::abs((likelihood_all.back() - likelihood_all[likelihood_all.size() - 2]) / likelihood_all[likelihood_all.size() - 2]);
      if (llk_diff < tol) {
        converge = true;
        break;
      }
    }

    if (key == maxit) break;
  }

  return List::create(Named("beta") = beta,
                      Named("likelihood_all") = likelihood_all,
                      Named("j_star_all") = j_star_all,
                      Named("converge") = converge);
}





