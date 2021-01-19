#ifndef PG_HPP
#define PG_HPP

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::colvec rpg(arma::colvec shape, arma::colvec scale);

#endif
