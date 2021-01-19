#include "rndpp_mvnormal.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]

arma::mat rndpp_mvnormal(int n, const arma::vec &mean, const arma::mat &sigma){
  int dimension = mean.n_elem;
  arma::mat outmat = arma::zeros(dimension, n);
  arma::mat cholsigma = arma::chol(sigma, "lower");
  arma::mat xtemp = (arma::randn(n, dimension)).t();
  arma::mat term2 = cholsigma * xtemp;
  for(int i=0; i<n; i++){
    outmat.col(i) = mean + term2.col(i);
  }
  return outmat.t();
}
