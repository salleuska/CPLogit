#ifndef DISTFROM_HPP
#define DISTFROM_HPP

// Using code from fxt library
// CREDITS TO GIVE
#include  "setpart.h"
#include  "fxtio.h"
#include  "fxttypes.h"
#include  "nextarg.h"
#include <math.h>

#include "vi_dist.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

List dist_from(arma::vec clust,bool return_partitions );

#endif
