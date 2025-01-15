#ifndef _F_C_H
#define _F_C_H


#include "f_mmult.h"
#include "f_integral.h"
#include "f_cgb.h"

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List f_c(NumericVector cgb, NumericVector Queue, int SNT, NumericVector OT_label, NumericVector B, double c0, double sigma_c0, NumericVector gb, NumericVector DELTA, NumericVector z, NumericVector OT, NumericVector Fail_label, NumericVector J, NumericVector at_c, NumericVector beta, NumericVector lambda, int N, int NY, int NZ, NumericVector c, int K, double c_c, int NP, int IF_QUADRATIC, NumericVector kns){
  NumericVector c_out = c;
  NumericMatrix beta1(1,NZ);
  NumericMatrix z1(NZ,1);
  
  for(int i=0; i<NZ; ++i){
    beta1(0,i) = beta[i];
  }
  
  for(int k=0; k<NY; ++k){
    //------------------------------ sample ------------------------------------//
    double c_star = rnorm(1,c[k],c_c)[0];
    //------------------------------ ratio ------------------------------------//
    double ratio_c = pow(c_star-c0,2)/(-2) - pow(c[k]-c0,2)/(-2 * pow(sigma_c0,2));
    
    for(int i=0; i<N; ++i){
      
      NumericVector u1(NP);
      NumericVector u1_star(NP);
      
      for(int p=0; p<NP; ++p){
        u1[p] = cgb[i+p*N];
        for(int s=0; s<NY; ++s){
          if(s == k){
            u1_star[p] += c_star*gb[i+(s*NP+p)*N];
          }else{
            u1_star[p] += c[s]*gb[i+(s*NP+p)*N];
          }
        }
      }
      
      if(DELTA[i]==1){
        for(int p=0; p<NP; ++p){
          ratio_c += B[p*SNT+Queue[i]+OT_label[i]]*(u1_star[p] - u1[p]);
        }
      }
      
      for(int j=0; j<NZ; ++j){
        z1(j,0) = z[j*N+i];
      }
      
      ratio_c += (f_integral(u1, J, lambda, Fail_label[i], K, OT[i], 1, 0, Fail_label[i], NP, IF_QUADRATIC, kns)-f_integral(u1_star, J, lambda, Fail_label[i], K, OT[i], 1, 0, Fail_label[i], NP, IF_QUADRATIC, kns))*exp(f_mmult(beta1,z1)(0,0));
      
    }
    //----------------------------- judge -----------------------------------//
    ratio_c = exp(ratio_c);
    
    double randnum = runif(1)[0];
    if(randnum<ratio_c){
      c_out[k] = c_star;
      at_c[k] += 1;
    }
    cgb = f_cgb(c_out, gb, N, NY, NP);
    
  }
  
  return List::create(Rcpp::Named("c") = c_out,
                      Rcpp::Named("at_c") = at_c);
}


#endif
