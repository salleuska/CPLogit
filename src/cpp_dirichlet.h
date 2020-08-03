#ifndef DIRICHLET_DIST_HPP
#define DIRICHLET_DIST_HPP

#include <Rcpp.h>
#include "shared.h"

// using std::pow;
// using std::sqrt;
// using std::abs;
// using std::exp;
// using std::log;
// using std::floor;
// using std::ceil;

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

NumericVector cpp_ddirichlet(
    const NumericMatrix& x,
    const NumericMatrix& alpha,
    const bool& log_prob
);

NumericMatrix cpp_rdirichlet(
    const int& n,
    const NumericMatrix& alpha
);

#endif
