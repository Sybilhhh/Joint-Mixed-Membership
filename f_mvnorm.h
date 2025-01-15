#ifndef _F_MVNORM_H
#define _F_MVNORM_H


 
using namespace arma; 
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]



NumericMatrix f_mvnorm(int n, NumericVector mu, NumericMatrix sigma) {
  arma::mat sigma_temp= Rcpp::as<arma::mat>(sigma);
  arma::vec mu_temp = Rcpp::as<arma::vec>(mu);
  int ncols = sigma_temp.n_cols;
  mat Y = randn(n, ncols);
  return(wrap((repmat(mu_temp, 1, n).t()+Y * chol(sigma_temp)).t()));
}

#endif
