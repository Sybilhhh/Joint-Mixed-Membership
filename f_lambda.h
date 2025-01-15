#ifndef _F_LAMBDA_H
#define _F_LAMBDA_H

#include "f_mmult.h"
#include "f_integral.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List f_lambda(NumericVector J, NumericVector lambda, NumericVector z, NumericVector OT, NumericVector DELTA, NumericVector cgb, NumericVector Fail_label, NumericVector beta, int N, int NY, int K, int NZ, double alpha1, double alpha2, NumericVector c, int NP, int IF_QUADRATIC, NumericVector kns) {
  NumericVector lambda_out = lambda;
  int L = J.size()-1; // number of intervals for estimating baseline hazard rate
  NumericVector lambda_null(L);
  NumericMatrix beta1(1,NZ);
  NumericVector sum1(L);
  NumericVector sum2(L);
  for(int i=0; i<L; ++i){
    lambda_null[i] = 1;
  }
  for(int i=0; i<NZ; ++i){
    beta1(0,i) = beta[i];
  }
  
  for(int i=0; i<L; ++i){
    for(int j=0; j<N; ++j){
      NumericMatrix z1(NZ,1);

      NumericVector u1(NP);
      for(int p=0; p<NP; ++p){
        u1[p] = cgb[j+p*N];
      }
      for(int k=0; k<NZ; ++k){
        z1(k,0) = z[k*N+j];
      }
      
      if((Fail_label[j]==(i+1))&(DELTA[j]==1)){
        sum1[i] += 1;
      }
      if(Fail_label[j]>=(i+1)){
        sum2[i] += f_integral(u1, J, lambda_null, Fail_label[j], K, OT[j], 1, i, i+1, NP, IF_QUADRATIC, kns)*exp(f_mmult(beta1,z1)(0,0));
      }
            
    }
    lambda_out[i] = rgamma(1,alpha1+sum1[i], 1/(alpha2+sum2[i]))[0];
  }
  return List::create(Rcpp::Named("lambda") = lambda_out,
                      Rcpp::Named("sum1") = sum1,
                      Rcpp::Named("sum2") = sum2);  
}
#endif
