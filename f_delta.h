#ifndef _F_DELTA_H
#define _F_DELTA_H

#include "f_mmult.h"

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List f_delta(double delta_0, NumericVector xi, NumericVector g, int N, int M, double a0_delta, double b0_delta, double c_delta, double at_delta){
  double delta_0_out = delta_0;
  NumericVector xi_out = xi;
  NumericVector delta(M);
  NumericVector delta_star(M);
  NumericMatrix delta_diff(1,M);
  NumericMatrix g1(M,1);
  
  // sample new delta_star
  for(int m=0; m<M; ++m){
    delta[m] = delta_0 * xi[m];
    delta_star[m] = rlnorm(1, log(delta[m]), c_delta)[0];
  }
  double delta_star_0 = sum(delta_star);

  for(int m=0; m<M; ++m){
    delta_diff(0,m) = delta_star[m] - delta[m];
    for(int i=0; i<N; ++i){
      g1(m,0) += log(g[m*N+i]);
    }
  }
  
  // compute ratio
  double ratio_1 = b0_delta * (delta_0 - delta_star_0);
  double ratio_2 = sum(log(delta_star) - log(delta));
  double ratio_3 = N * (lgamma(delta_star_0) - lgamma(delta_0));
  double ratio_4 = 0;
  for(int m=0; m<M; ++m){
    ratio_4 += (lgamma(delta[m]) - lgamma(delta_star[m]));
  }
  ratio_4 = N *ratio_4;
  double ratio_5 = f_mmult(delta_diff,g1)(0,0);
  double ratio_6 = (a0_delta - 1) * (log(delta_star_0) - log(delta_0));
  double ratio_delta = ratio_1 + ratio_2 + ratio_3 + ratio_4 + ratio_5 + ratio_6;

  // judge
  ratio_delta = exp(ratio_delta);
  
  double randnum = runif(1)[0];
  if(randnum<ratio_delta){
    delta_0_out = delta_star_0;
    xi_out = delta_star / delta_star_0;
    at_delta += 1;
  }
  
  return List::create(Rcpp::Named("delta_0") = delta_0_out,
                      Rcpp::Named("xi") = xi_out,
                      Rcpp::Named("at_delta") = at_delta,
                      Rcpp::Named("ratio_delta") = ratio_delta
  );
}
#endif
