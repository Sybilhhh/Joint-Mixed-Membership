#ifndef _F_SWITCHING_H
#define _F_SWITCHING_H

#include <f_judge.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List f_switching(NumericVector b, NumericVector xi, NumericVector g, int N, int M, int NY, int NP) {
  NumericVector b_out = clone(b);
  NumericVector xi_out = clone(xi);
  NumericVector g_out = clone(g);
  
  NumericVector mu(M);
  for (int m = 0; m < M; ++m){
    mu[m] = b[m*NY*NP+0*NP+0];
  }
  
  List all = f_judge(mu, M);
  bool judge = all[0];
  IntegerVector norder = all[1];
  
  if(judge){
    for(int m=0; m<M; ++m){
      int a = norder[m]-1;
      for(int k=0; k<NY; ++k){
        for (int p=0; p<NP; ++p){
          b_out[m*NY*NP+k*NP+p] = b[a*NY*NP+k*NP+p];
        }
      }
      xi_out[m] = xi[a];
      for(int i=0; i<N; ++i){
        g_out[m*N+i] = g[a*N+i];
      }
    }
    
  }
  
  return List::create(Rcpp::Named("b") = b_out,
                      Rcpp::Named("xi") = xi_out,
                      Rcpp::Named("g") = g_out
  );
}

#endif
