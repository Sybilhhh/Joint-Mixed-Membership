#ifndef _F_CGB_H
#define _F_CGB_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector f_cgb(NumericVector c, NumericVector gb, int N, int NY, int NP) {
  NumericVector cgb_out(N*NP);

  for(int i=0; i<N; ++i){
    for(int p=0; p<NP; ++p){
      for(int k=0; k<NY; ++k){
        cgb_out[i+p*N] += c[k]*gb[i+(k*NP+p)*N];
      }
    }
  }
  return cgb_out;
}
#endif
