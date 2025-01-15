#ifndef _F_INV_H
#define _F_INV_H

#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

NumericMatrix f_inv(NumericMatrix tmm){
  
  const Eigen::Map<Eigen::MatrixXd> ttm(as<Eigen::Map<Eigen::MatrixXd> >(tmm));
  
  Eigen::MatrixXd inv = ttm.inverse();
  return(wrap(inv));
}

#endif
