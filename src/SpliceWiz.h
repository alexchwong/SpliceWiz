#ifndef CODE_SPLICEWIZ
#define CODE_SPLICEWIZ

// Declare SPLICEWIZ in Makevars
// If not declared, this will produce an executable that does not depend on Rcpp
#ifdef SPLICEWIZ
  #include <Rcpp.h>
  using namespace Rcpp;
  
  // [[Rcpp::depends(RcppProgress)]]
  #include <progress.hpp>
  
  #define cout Rcpp::Rcout
#else
  #include <iostream>
#endif

#endif