#ifndef _F_INTEGRAL_H
#define _F_INTEGRAL_H


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double f_integral(NumericVector u, NumericVector J, NumericVector lambda, int Fail_label, int K, double OT, double c, int start, int end, int NP, int IF_QUADRATIC, NumericVector kns) {
  double integral =0;
  // Rcout << "The value of second_level : " << second_level << "\n";
  for(int i=start; i<end; ++i){
    
    if(i==(Fail_label-1)){
      for(int j=0; j<K; ++j){
        double t = J[i]+(OT-J[i])*j/K;
        double temp = u[0];
        if(NP>1){
          temp += u[1]*t;
          if(NP>2){
            if(IF_QUADRATIC==1){
              temp += u[2]*pow(t,2);
            }else{
              for(int p=2; p<NP; ++p){
                temp += u[p]*(-1)*((std::abs(pow(t-kns[NP-1],3))+pow(t-kns[NP-1],3)-(std::abs(pow(t-kns[p-2],3))+pow(t-kns[p-2],3)))/(2*(kns[NP-1]-kns[p-2]))-(std::abs(pow(t-kns[NP-1],3))+pow(t-kns[NP-1],3)-(std::abs(pow(t-kns[NP-2],3))+pow(t-kns[NP-2],3)))/(2*(kns[NP-1]-kns[NP-2])));
              }
            }
          }
        }
        // Rcout << "The value of temp : " << temp << "\n";
        integral += lambda[i]*exp(c*temp)*(OT-J[i])/K;
      }
    }else{
      for(int j=0; j<K; ++j){
        double t = J[i]+(J[i+1]-J[i])*j/K;
        double temp = u[0];
        if(NP>1){
          temp += u[1]*t;
          if(NP>2){
            if(IF_QUADRATIC==1){
              temp += u[2]*pow(t,2);
            }else{
              for(int p=2; p<NP; ++p){
                temp += u[p]*(-1)*((std::abs(pow(t-kns[NP-1],3))+pow(t-kns[NP-1],3)-(std::abs(pow(t-kns[p-2],3))+pow(t-kns[p-2],3)))/(2*(kns[NP-1]-kns[p-2]))-(std::abs(pow(t-kns[NP-1],3))+pow(t-kns[NP-1],3)-(std::abs(pow(t-kns[NP-2],3))+pow(t-kns[NP-2],3)))/(2*(kns[NP-1]-kns[NP-2])));
              }
            }
          }
        }
        
        // Rcout << "The value of temp : " << temp << "\n";
        
        integral += lambda[i]*exp(c*temp)*(J[i+1]-J[i])/K;
      }
    }
    
  }
  
  return integral;  
}

#endif
