#ifndef _F_GB_H
#define _F_GB_H

#include "f_mmult.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector f_gb(NumericVector b, NumericVector g, int N, int M, int NY, int NP) {
  NumericVector gb_out(N*NY*NP);
  NumericMatrix g1(1,M);
  NumericMatrix b1(M,1);
  NumericVector cgb_out(N*NY*NP);
  
  for(int j=0; j<(NY*NP); ++j){
    for(int m=0; m<M; ++m){
      b1(m,0) = b[m*NY*NP+j];
    }
    for(int i=0; i<N; ++i){
      for(int m=0; m<M; ++m){
        g1(0,m) = g[m*N+i];
      }
      gb_out[i+j*N] = f_mmult(g1,b1)(0,0);
    }
  }
  return gb_out;
}
#endif
