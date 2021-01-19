// Function to sample partitions uniformly
#ifndef SAMPLEPARTITION_HPP
#define SAMPLEPARTITION_HPP

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

List samplePartition();

#endif
