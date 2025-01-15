#ifndef _F_BETA_H
#define _F_BETA_H

#include "f_inv.h"
#include "f_mmult.h" 
#include "f_madd.h"
#include "f_mvnorm.h"
#include "f_integral.h"

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List f_beta(NumericVector beta0, NumericMatrix Sigma_beta0, NumericVector cgb, NumericVector DELTA, NumericVector z, NumericVector OT, NumericVector Fail_label, NumericVector J, int at_beta, NumericVector beta, NumericVector lambda, int N, int NY, int NZ, NumericVector c, int K, double c_beta, int NP, int IF_QUADRATIC, NumericVector kns){
  NumericVector beta_out = beta;
  NumericMatrix sigma_beta = f_inv(Sigma_beta0);
  NumericMatrix z1(NZ,1);
  NumericMatrix z2(1,NZ);
  NumericMatrix beta_star1(1,NZ);
  NumericMatrix beta_star2(NZ,1);
  NumericMatrix beta01(1,NZ);
  NumericMatrix beta02(NZ,1);
  NumericMatrix beta_old1(1,NZ);
  NumericMatrix beta_old2(NZ,1);
  NumericVector integral(N);
  NumericVector u1(NP);


  for(int i=0; i<N; ++i){
    for(int j=0; j<NZ; ++j){
      z1(j,0) = z[j*N+i];
      z2(0,j) = z[j*N+i];
    }
    
    for (int p = 0; p < NP; ++p){
      u1[p] = cgb[i+p*N];
    }
    integral[i] = f_integral(u1, J, lambda, Fail_label[i], K, OT[i], 1, 0, Fail_label[i], NP, IF_QUADRATIC, kns);
    // Rcout << "The value of integral : " << integral[i] << "\n";
    sigma_beta = f_madd(sigma_beta,f_mmult(z1,z2)*integral[i]);
  }

  sigma_beta = f_inv(sigma_beta);

  beta_star2 = f_mvnorm(1,beta,sigma_beta*c_beta);

  for(int i=0; i<NZ; ++i){
    beta_old1(0,i) = beta[i];
    beta_old2(i,0) = beta[i];
    beta01(0,i) = beta0[i];
    beta02(i,0) = beta0[i];
    beta_star1(0,i) = beta_star2(i,0);
  }

  double ratio_beta = f_mmult(f_mmult(f_madd(beta_star1,beta01*(-1)),f_inv(Sigma_beta0)),f_madd(beta_star2,beta02*(-1)))(0,0)/(-2)-f_mmult(f_mmult(f_madd(beta_old1,beta01*(-1)),f_inv(Sigma_beta0)),f_madd(beta_old2,beta02*(-1)))(0,0)/(-2);

  for(int i=0; i<N; ++i){
    for(int j=0; j<NZ; ++j){
      z1(j,0) = z[j*N+i];
    }
    ratio_beta += f_mmult(f_madd(beta_star1,beta_old1*(-1)),z1)(0,0)*DELTA[i]-integral[i]*(exp(f_mmult(beta_star1,z1)(0,0))-exp(f_mmult(beta_old1,z1)(0,0)));
  }

  ratio_beta = exp(ratio_beta);

  double randnum = runif(1)[0];
  if(randnum<ratio_beta){
    for(int i=0; i<NZ; ++i){
      beta_out[i] = beta_star2(i,0);
    }
    at_beta += 1;
  }
  

  
  return List::create(Rcpp::Named("beta") = beta_out,
                      Rcpp::Named("at_beta") = at_beta,
                      Rcpp::Named("ratio_beta") = ratio_beta,
                      Rcpp::Named("beta_star2") = beta_star2,
                      Rcpp::Named("beta_star1") = beta_star1,
                      Rcpp::Named("beta_old1") = beta_old1,
                      Rcpp::Named("Sigma_beta0") = Sigma_beta0,
                      Rcpp::Named("integral") = integral);
}
#endif
