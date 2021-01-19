#ifndef UTILS_HPP
#define UTILS_HPP

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// Function which normalize probabilities in log scale
arma::vec normalize_probs(arma::vec x);

arma::uvec relabel(arma::uvec current);

#endif
