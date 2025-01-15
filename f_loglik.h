#ifndef _F_LOGLIK_H
#define _F_LOGLIK_H

#include "f_mmult.h" 
#include "f_integral.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double f_loglik(NumericVector cgb, NumericVector OT_label, NumericVector B, NumericVector J, NumericVector Fail_label, NumericVector z, NumericVector OT, NumericVector DELTA, NumericVector Queue, NumericVector y, NumericVector NT, double pi, int K, int M, int NY, int SNT, int N, int NZ, NumericVector lambda, NumericVector beta, NumericVector c, NumericVector g, NumericVector b, NumericVector psi_eps, NumericVector gb, int NP, int IF_QUADRATIC, NumericVector kns) {
  
  double dic_1 = 0;
  double dic_2 = 0;
  double dic_3 = 0;
  double dic_4 = 0;
  double dic_5 = 0;
  
  NumericMatrix beta1(1,NZ);
  for(int i=0; i<NZ; ++i){
    beta1(0,i) = beta[i];
  }
  
  
  for(int i=0; i<N; ++i){
    for(int k=0; k<NY; ++k){
      // --------------- MMM ------------------------
      int t_temp = NT[i];
      for(int j=0; j<t_temp; ++j){
        double Y_n = y[Queue[i]+j+k*SNT];
        for (int p = 0; p < NP; ++p){
          Y_n -= gb[(k*NP+p)*N+i]*B[p*SNT+Queue[i]+j];
        }
        dic_1 -= pow(Y_n,2)/(2*psi_eps[k]);
        dic_1 -= 0.5* log(2*pi*psi_eps[k]);
      }
    }
    // -------------- COX --------------------------
    int L = J.size()-1;
    for(int j=0; j<L; ++j){
      if((Fail_label[i]==(j+1))&(DELTA[i]==1)){
        dic_2 += log(lambda[j]);
      }
    }
    // -------------- 3 ----------------------------
    NumericMatrix z1(NZ,1);
    for(int k=0; k<NZ; ++k){
      z1(k,0) = z[k*N+i];
    }
    dic_3 += DELTA[i] * f_mmult(beta1,z1)(0,0);
    // -------------- 4 ----------------------------
    NumericVector u1(NP);
    for(int p=0; p<NP; ++p){
      u1[p] = cgb[i+p*N];
      dic_4 += DELTA[i]*B[p*SNT+Queue[i]+OT_label[i]]*u1[p];
    }
    dic_5 += exp(f_mmult(beta1,z1)(0,0))*f_integral(u1, J, lambda, Fail_label[i], K, OT[i], 1, 0, Fail_label[i], NP, IF_QUADRATIC, kns);
  }
  double dic = dic_1 + dic_2 + dic_3 + dic_4 + dic_5;
  
  return dic;
  
}

#endif
