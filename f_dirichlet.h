#ifndef _F_DIRICHLET_H
#define _F_DIRICHLET_H

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]

NumericVector f_dirichlet(NumericVector a, int K){
  NumericVector y_out(K);
  NumericVector x_temp(K);
  double x_sum = 0;
  
  for(int k = 0; k<K; ++k){
    x_temp[k] = rgamma(1, a[k], 1)[0];
    x_sum += x_temp[k];
  }
  
  for(int k =0; k<K; ++k){
    y_out[k] = x_temp[k]/x_sum;
  }
  
  return y_out;
}

#endif
