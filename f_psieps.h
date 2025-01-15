#ifndef _F_PSIEPS_H
#define _F_PSIEPS_H

#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
NumericVector f_psieps(NumericVector psi_eps,NumericVector y, NumericVector gb, NumericVector B, NumericVector Queue, NumericVector NT, int N, int SNT, int NY, double alpha_eps0, double beta_eps0, int NP) {
  NumericVector psieps_out = psi_eps;
  double alpha = alpha_eps0 + 0.5 * SNT;
  
  for(int k=0; k<NY; ++k){
    double sum = 0;
    
    for(int i=0; i<N; ++i){
      int t_temp = NT[i];
      
      for(int j=0; j<t_temp; ++j){
        double h_temp = 0;
        for(int p=0;p<NP; ++p){
          h_temp +=  gb[(k*NP+p)*N+i]*B[p*SNT+Queue[i]+j]; 
        } 
        double temp = y[Queue[i]+j+k*SNT] - h_temp;
        sum += pow(temp,2);
      }
    }
    double beta = beta_eps0 + 0.5 *sum;
    psieps_out[k] = 1/(rgamma(1,alpha,(1/beta))[0]);
    // Rprintf("alpha: %f; beta: %f\n",alpha, beta);
  }
  
  return psieps_out;
}


#endif
