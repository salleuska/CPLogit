#ifndef MARGINALVI_HPP
#define MARGINALVI_HPP

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec marginal_vi(arma::uvec clustering, arma::uvec C_0, int H, int obs);

#endif
