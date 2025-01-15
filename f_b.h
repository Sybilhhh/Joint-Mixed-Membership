#ifndef _F_B_H
#define _F_B_H

#include "f_inv.h"
#include "f_mmult.h" 
#include "f_madd.h"
#include "f_integral.h"
#include "f_gb.h"
#include "f_cgb.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List f_b(NumericVector cgb, NumericVector gb, NumericVector b, NumericVector y, NumericVector B, NumericVector g, NumericVector psi_eps, NumericVector OT, NumericVector DELTA, NumericVector Fail_label, NumericVector J, NumericVector lambda, NumericVector beta, NumericVector z, NumericVector Queue, NumericVector c, NumericVector NT, int N, int M, int NY, int K, int NZ, int SNT, double c_b, NumericVector at_b, NumericVector b0, NumericMatrix Sigma_b0, int NP, int IF_QUADRATIC, NumericVector kns, NumericVector OT_label){
  NumericVector b_out = b;
  NumericMatrix beta1(1,NZ);
  NumericMatrix z1(NZ,1);
  NumericMatrix bmk_star1(1,NP);
  NumericMatrix bmk_star2(NP,1);
  NumericMatrix bmk_old1(1,NP);
  NumericMatrix bmk_old2(NP,1);
  NumericMatrix bmk01(1,NP);
  NumericMatrix bmk02(NP,1);
  
  for(int i=0; i<NZ; ++i){
    beta1(0,i) = beta[i];
  }
  
  for(int m=0; m<M; ++m){
    for(int k=0; k<NY; ++k){
      //------------------------------ sample ------------------------------------//
      for (int p = 0; p < NP; ++p){
        double temp_Sigma_b = 0;
        for(int i=0; i<N; ++i){
          int t_temp = NT[i];
          for(int j=0; j<t_temp; ++j){
            temp_Sigma_b += pow((g[m*N+i]*B[p*SNT+Queue[i]+j]),2);
          }
        }
        double sigma_bmk = 1/(1/Sigma_b0(p,p)+temp_Sigma_b/psi_eps[k]);
        bmk_star2(p,0) =  rnorm(1,b[m*NY*NP+k*NP+p],sqrt(c_b * sigma_bmk))[0];
        bmk_old1(0,p) = b[m*NY*NP+k*NP+p];
        bmk_old2(p,0) = b[m*NY*NP+k*NP+p];
        bmk01(0,p) = b0[p];
        bmk02(p,0) = b0[p];
        bmk_star1(0,p) = bmk_star2(p,0);
      }      
      //------------------------------ ratio ------------------------------------//
      double ratio_1 = f_mmult(f_mmult(f_madd(bmk_old1,bmk01*(-1)),f_inv(Sigma_b0)),f_madd(bmk_old2,bmk02*(-1)))(0,0)/2-f_mmult(f_mmult(f_madd(bmk_star1,bmk01*(-1)),f_inv(Sigma_b0)),f_madd(bmk_star2,bmk02*(-1)))(0,0)/2;
      double ratio_2 = 0;
      double ratio_3 = 0;
      double ratio_4 = 0;
      
      for(int i=0; i<N; ++i){
        NumericVector gb_star(NP);
        for (int j = 0; j < NP; ++j){
          for (int s = 0; s < M; ++s){
            if (s == m){
              gb_star[j] += bmk_star2(j,0) * g[s*N+i];
            }else{
              gb_star[j] += b[s*NY*NP+k*NP+j] * g[s*N+i];
            }
          }
        }
        int t_temp = NT[i];
        for(int j=0; j<t_temp; ++j){
          double Y_n = y[Queue[i]+j+k*SNT];
          double Y_star = y[Queue[i]+j+k*SNT];
          for (int s = 0; s < NP; ++s){
            Y_n -= gb[(k*NP+s)*N+i]*B[s*SNT+Queue[i]+j];
            Y_star -= gb_star[s]*B[s*SNT+Queue[i]+j];
          }
          ratio_2 += (pow(Y_n,2) - pow(Y_star,2))/(2*psi_eps[k]);
        }
        //------------------------------ cox ------------------------------------//
        for (int p = 0; p < NP; ++p){
          ratio_3 += DELTA[i] * B[p*SNT+Queue[i]+OT_label[i]] * c[k] * g[m*N+i]*(bmk_star2(p,0) - b[m*NY*NP+k*NP+p]);
        }
        //------------------------------ ratio 4 ------------------------------------//
        NumericVector u1(NP);
        NumericVector u1_star(NP);

        for (int j = 0; j < NP; ++j){
          u1[j] = cgb[i+j*N];
        }

        for (int j = 0; j < NP; ++j){
          for (int t = 0; t < NY; ++t){
            for (int s = 0; s < M; ++s){
              if (s == m && t == k){
                u1_star[j] += c[t]* bmk_star2(j,0) * g[s*N+i];
              }else{
                u1_star[j] += c[t]* b[s*NY*NP+t*NP+j] * g[s*N+i];
              }
            }
          }
        }
      
        for(int j=0; j<NZ; ++j){
          z1(j,0) = z[j*N+i];
        }
        ratio_4 += exp(f_mmult(beta1,z1)(0,0)) *(f_integral(u1, J, lambda, Fail_label[i], K, OT[i], 1, 0, Fail_label[i], NP, IF_QUADRATIC, kns) - f_integral(u1_star, J, lambda, Fail_label[i], K, OT[i], 1, 0, Fail_label[i], NP, IF_QUADRATIC, kns));
      }
        
      //----------------------------- judge -----------------------------------//
      double log_ratio = ratio_1 + ratio_2 + ratio_3 + ratio_4;
      double ratio_bi = exp(log_ratio);
      double randnum = runif(1)[0];
      if(randnum<ratio_bi){
        for (int p = 0; p < NP; ++p){
          b_out[m*NY*NP+k*NP+p] = bmk_star2(p,0);
        }
        at_b[m*NY+k] += 1;
      }
      gb = f_gb(b_out, g, N, M, NY, NP);
      cgb = f_cgb(c, gb, N, NY, NP);
      //Rprintf("Sigma_b: %f;bmk_star: %f\n",sigma_bmk, bmk_star);
      //Rprintf("ratio1: %f; ratio2: %f; ratio3: %f; ratio4: %f\n",ratio_1,ratio_2,ratio_3,ratio_4);  
    }
  }
  return List::create(Rcpp::Named("b") = b_out,
                      Rcpp::Named("at_b") = at_b
                      //,Rcpp::Named("gb") = gb,
                      //Rcpp::Named("cgb") = cgb
  );
}

#endif
