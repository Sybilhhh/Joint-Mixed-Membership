#ifndef _F_JUDGE_H
#define _F_JUDGE_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List f_judge(NumericVector mu, int M) {
  
  if (is_true(any(duplicated(mu)))) {
    Rf_warning("There are duplicates in 'mu'; order not guaranteed to match that of R's base::order");
  }
  NumericVector sorted = clone(mu).sort();
  IntegerVector norder = match(sorted, mu);
  
  bool judge = FALSE;
  for(int m=0; m<M; ++m){
    if(norder[m] != m+1){
      judge = TRUE;
    }
  }
  return List::create(Rcpp::Named("judge") = judge,
                      Rcpp::Named("norder") = norder);
}

#endif
