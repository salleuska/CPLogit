// Polya gamma random generator
#include "RcppArmadillo.h"
#include <R_ext/Utils.h>
#include <iostream>
#include "PolyaGamma.h"
#include "PG.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::colvec rpg(arma::colvec shape, arma::colvec scale) {
// C++-only interface to PolyaGamma class
// draws random PG variates from arma::vectors of n's and psi's
  RNG r;
  PolyaGamma pg;
#ifdef USE_R
  GetRNGstate();
#endif
  int d = shape.n_elem;
  arma::colvec result(d);
  for(int i=0; i<d; i++) {
    result[i] = pg.draw(shape(i), scale(i), r);
  }
#ifdef USE_R
  PutRNGstate();
#endif
  return result;
}
