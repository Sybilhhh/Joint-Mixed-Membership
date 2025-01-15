#ifndef _F_GI_H
#define _F_GI_H

#include "f_mmult.h"
#include "f_dirichlet.h"
#include "f_integral.h"

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List f_gi(NumericVector cgb, double delta_0, NumericVector xi, NumericVector gb, NumericVector b, NumericVector y, NumericVector B, NumericVector g, NumericVector psi_eps, NumericVector OT, NumericVector DELTA, NumericVector Fail_label, NumericVector J, NumericVector lambda, NumericVector beta, NumericVector z, NumericVector Queue, NumericVector c, NumericVector NT, int N, int M, int NY, int K, int NZ, int SNT, double d0, double c_g, NumericVector at_g, int NP, int IF_QUADRATIC, NumericVector kns, NumericVector OT_label){
  NumericVector g_out = g;
  NumericMatrix bB(M,1);
  NumericMatrix g1(1,M);
  NumericMatrix gi_star(1,M);
  NumericMatrix beta1(1,NZ);
  NumericVector delta(M);
  NumericMatrix delta1(M,1);
  NumericMatrix gi_diff(1,M);
  NumericMatrix z1(NZ,1);
  
  for(int m=0; m<M; ++m){
    delta[m] = c_g * d0;
    delta1(m,0) = delta_0 * xi[m] - 1;
  }
  for(int i=0; i<NZ; ++i){
    beta1(0,i) = beta[i];
  }
  
  
  for(int i=0; i<N; ++i){
    //------------------------------ sample ------------------------------------//
    NumericVector g_star = f_dirichlet(delta,M);
    
    //------------------------------ ratio 1------------------------------------//
    for(int j=0; j<M; ++j){
      if(g_star[j] == 0){
        g_star[j] = 1e-323;
      }
      if(g[j*N+i] == 0){
        g[j*N+i] = 1e-323;
      }
      g1(0,j) = g[j*N+i];
      gi_star(0,j) = g_star[j];
      gi_diff(0,j) = log(gi_star[j]) - log(g[j*N+i]);
    }
    double ratio_1 = f_mmult(gi_diff,delta1)(0,0);
    //------------------------------ ratio 2-----------------------------------//
    double ratio_2 = 0;
    int t_temp = NT[i];
    for(int k=0; k<NY; ++k){
      for(int j=0; j<t_temp; ++j){
        /*
        for(int m =0; m<M;++m){
          x(m,0) = mu[m+k*M] + b[k*M+m]* B[Queue[i]+j];
        }
        double Y_n = y[Queue[i]+j+k*SNT] - gmu[i+k*N]- gb[i+k*N] *B[Queue[i]+j];
        double Y_star = y[Queue[i]+j+k*SNT] - f_mmult(gi_star,x)(0,0);
        */
        double Y_n = y[Queue[i]+j+k*SNT];
        double Y_star = y[Queue[i]+j+k*SNT];
        for (int p = 0; p < NP; ++p){
          for(int m =0; m<M;++m){
            bB(m,0) = b[m*NY*NP+k*NP+p]*B[p*SNT+Queue[i]+j];
          }
          Y_n -= gb[(k*NP+p)*N+i]*B[p*SNT+Queue[i]+j];
          Y_star -= f_mmult(gi_star,bB)(0,0);
        }
        
        ratio_2 += (pow(Y_n,2) - pow(Y_star,2))/(2 * psi_eps[k]);
      }
    }
    //------------------------------ ratio 3-----------------------------------//
    /*
    NumericMatrix mu1(M,1);
    NumericMatrix b1(M,1);
    for(int k=0; k<NY; ++k){
      for(int m=0; m<M; ++m){
        mu1(m,0) += c[k]*mu[k*M+m];
        b1(m,0) += c[k]*b[k*M+m];
      }
    }
    NumericVector u1(2);
    NumericVector u1_star(2);
    u1[0] = f_mmult(g1,mu1)(0,0);
    u1_star[0] = f_mmult(g1,mu1)(0,0);
    u1[1] = f_mmult(g1,b1)(0,0);
    u1_star[1] = f_mmult(gi_star,b1)(0,0);
    double ratio_3 = DELTA[i] * OT[i]  * (u1_star[1] - u1[1]);
    ratio_3 += DELTA[i] * (u1_star[0] - u1[0]);
    */
    NumericMatrix b1(M,1);
    NumericVector u1(NP);
    NumericVector u1_star(NP);
    double ratio_3 = 0;

    for (int p = 0; p < NP; ++p){
      u1[p] = cgb[i+p*N];
    }

    for (int p = 0; p < NP; ++p){
      for (int k=0; k<NY; ++k){
        for (int m = 0; m < M; ++m){
          b1(m,0) = b[m*NY*NP+k*NP+p];
        }
        u1_star[p] += c[k]*f_mmult(gi_star,b1)(0,0);
      }
      ratio_3 += DELTA[i] * B[p*SNT+Queue[i]+OT_label[i]]  * (u1_star[p] - u1[p]);  
    }
    
    //------------------------------ ratio 4-----------------------------------//
    for(int k=0; k<NZ; ++k){
      z1(k,0) = z[k*N+i];
    }
    double ratio_4 = exp(f_mmult(beta1,z1)(0,0)) *(f_integral(u1, J, lambda, Fail_label[i], K, OT[i], 1, 0, Fail_label[i], NP, IF_QUADRATIC, kns) - f_integral(u1_star, J, lambda, Fail_label[i], K, OT[i], 1, 0, Fail_label[i], NP, IF_QUADRATIC, kns));
    //------------------------------ ratio 5-----------------------------------//
    
    double ratio_5 = 0;
    for (int m = 0; m < M; ++m){
      ratio_5 += (delta[m]-1)*(log(g[m*N+i]) - log(g_star[m]));
    }
    
    //----------------------------- judge -----------------------------------//
    double log_ratio = ratio_1 + ratio_2 + ratio_3 + ratio_4 + ratio_5;
    // double log_ratio = ratio_1 + ratio_2 + ratio_3 + ratio_4;
    double ratio_gi = exp(log_ratio);
    double randnum = runif(1)[0];
    if(randnum<ratio_gi){
      for(int m = 0; m<M; ++m){
        g_out[m*N+i] = g_star[m];
      }
      at_g[i] += 1;
    }
  }
  
  return List::create(Rcpp::Named("g") = g_out,
                      Rcpp::Named("at_g") = at_g
  );
}
#endif
