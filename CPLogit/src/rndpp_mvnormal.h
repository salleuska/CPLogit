#ifndef RNDPP_MVNORMAL_HPP
#define RNDPP_MVNORMAL_HPP

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat rndpp_mvnormal(int n, const arma::vec &mean, const arma::mat &sigma);

#endif
