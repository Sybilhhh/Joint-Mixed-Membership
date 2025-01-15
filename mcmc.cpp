#include <RcppArmadillo.h>

#include "f_mmult.h"
#include "f_gb.h"
#include "f_delta.h"
#include "f_psieps.h"
#include "f_b.h"
#include "f_gi.h"
#include "f_switching.h"
#include "f_beta.h"
#include "f_c.h"
#include "f_lambda.h"
#include "f_loglik.h"
#include "f_cgb.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List mcmc(List observed_data, List para, List mh, List prior, int K, int iter, int rep, int T, int S, double pi, int IF_MMM, int IF_COX, int IF_AIC, int IF_QUADRATIC){

  // Cox
  double alpha1 = prior[0];
  double alpha2 = prior[1];
  NumericVector beta0 = prior[2];
  NumericMatrix Sigma_beta0 = prior[3];
  double c0 = prior[4];
  double sigma_c0 = prior[5];
  // MMM
  NumericVector b0 = prior[6];
  NumericMatrix Sigma_b0 = prior[7];
  double a0_delta = prior[8];
  double b0_delta = prior[9];
  double d0 = prior[10];
  double alpha_eps0 = prior[11];
  double beta_eps0 = prior[12];
  

  double c_beta = mh[0];
  double c_c = mh[1];
  double c_delta = mh[2];
  double c_b = mh[3];
  double c_g = mh[4];
  
  NumericVector lambda = para[0];
  NumericVector beta = para[1];
  NumericVector c = para[2];
  NumericVector b = para[3];
  NumericVector g = para[4];
  NumericVector psi_eps = para[5];
  double delta_0 = para[6];
  NumericVector xi = para[7];
  
  NumericVector OT_label = observed_data[0];
  NumericVector B = observed_data[1];
  NumericVector J =observed_data[2];
  NumericVector Fail_label = observed_data[3];
  NumericVector z = observed_data[4];
  NumericVector OT = observed_data[5];
  NumericVector DELTA = observed_data[6];
  NumericVector Queue = observed_data[7];
  NumericVector y = observed_data[8];
  NumericVector NT = observed_data[9];
  NumericVector kns = observed_data[10];
  
  // MMM
  int M = xi.size();
  int NY = psi_eps.size();
  int NP = b.size()/(M*NY);
  // Cox
  int SNT = B.size()/NP;
  int N = g.size()/M;
  int NZ = beta.size();
  int J_len = J.size();
  
  NumericVector gb(N*NY*NP);
  NumericVector cgb(N*NP);
  
  // Cox
  NumericVector lambda_str(S*(J_len-1));
  NumericVector beta_str(S*NZ);
  NumericVector c_str(S*NY);
  // MMM
  NumericVector b_str(S*NP*NY*M);
  NumericVector delta_0_str(S);
  NumericVector xi_str(S*M);
  NumericVector g_str(S*N*M);
  NumericVector psieps_str(S*NY);
  
  NumericVector AIC1_str(S);
  
  int at_beta = 0;
  NumericVector at_c(NY);
  int at_delta = 0;
  NumericVector at_b(M*NY);
  NumericVector at_g(N);
  
  gb = f_gb(b, g, N, M, NY, NP);
  cgb = f_cgb(c, gb, N, NY, NP);
  
  for(int t=0; t<T; ++t){
    if (IF_MMM==1){
      
      List b_all = f_b(cgb, gb, b, y, B, g, psi_eps, OT, DELTA, Fail_label, J, lambda, beta, z, Queue, c, NT, N, M, NY, K, NZ, SNT, c_b, at_b, b0, Sigma_b0, NP, IF_QUADRATIC, kns, OT_label);
      b = b_all[0];
      at_b= b_all[1];
      gb = f_gb(b, g, N, M, NY, NP);
      cgb = f_cgb(c, gb, N, NY, NP);
      
      List g_all = f_gi(cgb, delta_0, xi, gb, b, y, B, g, psi_eps, OT, DELTA, Fail_label, J, lambda, beta, z, Queue, c, NT, N, M, NY, K, NZ, SNT, d0, c_g, at_g,NP,IF_QUADRATIC,kns,OT_label);
      g = g_all[0];
      at_g = g_all[1];
      gb = f_gb(b, g, N, M, NY, NP);
      cgb = f_cgb(c, gb, N, NY, NP);
      
      psi_eps = f_psieps(psi_eps, y, gb, B, Queue, NT, N, SNT, NY, alpha_eps0, beta_eps0, NP);

      List delta_all = f_delta(delta_0, xi, g, N, M, a0_delta, b0_delta, c_delta, at_delta);
      delta_0 = delta_all[0];
      xi = delta_all[1];
      at_delta = delta_all[2];
      
      List per = f_switching(b, xi, g, N, M, NY,NP);
      b = per[0];
      xi = per[1];
      g = per[2];
      gb = f_gb(b, g, N, M, NY, NP);
      cgb = f_cgb(c, gb, N, NY, NP);
      
    }
    if(IF_COX==1){
      List beta_all = f_beta(beta0, Sigma_beta0, cgb, DELTA, z, OT, Fail_label, J, at_beta, beta, lambda, N, NY, NZ, c, K, c_beta, NP, IF_QUADRATIC, kns);
      beta = beta_all[0];
      at_beta = beta_all[1];
      
      List c_all = f_c(cgb, Queue, SNT, OT_label, B, c0, sigma_c0, gb, DELTA, z, OT, Fail_label, J, at_c, beta, lambda, N, NY, NZ, c, K, c_c, NP, IF_QUADRATIC, kns);
      c = c_all[0];
      at_c = c_all[1];
      cgb = f_cgb(c, gb, N, NY, NP);
      
      List lambda_all = f_lambda(J, lambda, z, OT, DELTA, cgb, Fail_label, beta, N, NY, K, NZ, alpha1, alpha2, c, NP, IF_QUADRATIC, kns);
      lambda = lambda_all[0];
    }
    
    if(t>(T-S-1)){
      int ind = t-(T-S);
      
      delta_0_str[ind] = delta_0;
      
      for(int i=0; i<M*NY*NP; ++i){
        b_str[ind*M*NY*NP+i] = b[i];
      }
      
      for(int i=0; i<M; ++i){
        xi_str[ind*M+i] = xi[i];
      }
      
      for(int i = 0; i<N*M; ++i){
        g_str[ind*N*M+i] = g[i];
      }
      
      for(int i=0; i<NY; ++i){
        psieps_str[ind*NY+i] = psi_eps[i];
        c_str[ind*NY+i] = c[i];
      }
      
      for(int i=0; i<(J_len-1); ++i){
        lambda_str[ind*(J_len-1)+i] = lambda[i];
      }
      for(int i=0; i<NZ; ++i){
        beta_str[ind*NZ+i] = beta[i];
      }
      
      // ------------- dic ------------------------------
      if(IF_AIC == 1){
        double temp_loglik = f_loglik(cgb, OT_label, B, J, Fail_label, z, OT, DELTA, Queue, y, NT, pi, K, M, NY, SNT, N, NZ, lambda, beta, c, g, b, psi_eps, gb, NP, IF_QUADRATIC, kns);
        AIC1_str[ind] = -2 * temp_loglik  + 2* (J_len-1+NZ+NY+M*NY*NP+M*N+NY+1+M);
      }
       
       
    }
    
    if((t+1)%100==0){
      Rprintf("\rRunning MCMC, %d/%d replication, %f%% completed...",iter,rep,(t+1)*100.0/T);
    }
    
  }
  
  return List::create(Rcpp::Named("lambda_str") = lambda_str,
                      Rcpp::Named("beta_str") = beta_str,
                      Rcpp::Named("c_str") = c_str,
                      Rcpp::Named("at_c") = at_c,
                      Rcpp::Named("at_beta") = at_beta,
                      Rcpp::Named("b_str") = b_str,
                      Rcpp::Named("at_b") = at_b,
                      Rcpp::Named("delta_0_str") = delta_0_str,
                      Rcpp::Named("xi_str") = xi_str,
                      Rcpp::Named("at_delta") = at_delta,
                      Rcpp::Named("g_str") = g_str,
                      Rcpp::Named("at_g") = at_g,
                      Rcpp::Named("psieps_str") = psieps_str
                      , Rcpp::Named("gb") = gb
                      , Rcpp::Named("AIC") = AIC1_str
                      , Rcpp::Named("cgb") = cgb
                      );
}