// Prior - Dirichlet Process
#ifndef PRIORDP_HPP
#define PRIORDP_HPP

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

class priorDP{
  // Prior distribution for mixture model weights - Dirichlet process
public:
  //---- Hyperparameters ----//
  // number of initial classes
  double H0;
  // upper bound for number of classes
  double H_up;
  // concentration parameter
  double alpha0;
  // Number of observations
  double n;
  // Log-robability of a new class
  double log_new_class;

public:
  //=========== CONSTRUCTORS =============//
  // Default constructor
  priorDP(double H, double upper_bound, double alpha, double n_obs)
  {
    H0 = H;
    n = n_obs;
    H_up = upper_bound;
    alpha0 = alpha;
    log_new_class = log(alpha0) - log(alpha0 + n - 1);
  }

  //=========== METHODS =============//
  // Log conditional distribution
  double log_conditional(double n_class);
  arma::vec log_conditional(arma::vec n_class);

};

//=========== METHODS =============//
// Conditional prior distribution in log scale
double priorDP::log_conditional(double n_class)
{
 return log(n_class) - log(alpha0 + n - 1);
}

// Conditional prior distribution in log scale - vector vers
arma::vec priorDP::log_conditional(arma::vec n_class)
{
 return log(n_class) - log(alpha0 + n - 1);
}

#endif
