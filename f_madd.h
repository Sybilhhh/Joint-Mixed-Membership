#ifndef _F_MADD_H
#define _F_MADD_H

#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

NumericMatrix f_madd(NumericMatrix tmm, NumericMatrix tm22){
  
  const Eigen::Map<Eigen::MatrixXd> ttm(as<Eigen::Map<Eigen::MatrixXd> >(tmm));
  const Eigen::Map<Eigen::MatrixXd> ttm2(as<Eigen::Map<Eigen::MatrixXd> >(tm22));
  
  Eigen::MatrixXd sum = ttm+ttm2;
  return(wrap(sum));
}

#endif 
