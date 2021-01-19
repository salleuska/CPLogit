#ifndef VIDIST_HPP
#define VIDIST_HPP

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// Function which normalize probabilities in log scale

double entropy(arma::uvec cl);

double vi_distC(arma::uvec cl1, arma::uvec cl2);

#endif
