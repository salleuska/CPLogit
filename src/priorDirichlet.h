// Pior - Dirichlet - Overfitted mixture

#ifndef PRIORDIRICHLET_HPP
#define PRIORDIRICHLET_HPP

#include <RcppArmadillo.h>
#include "includes.h"
#include "cpp_dirichlet.h"

// [[Rcpp::depends(RcppArmadillo)]]
 using namespace Rcpp;

class priorDirichlet{
  // Prior distribution for mixture model weights - symmetric Dirichlet
  // nu = mixture weights ~ Dirichlet(alpha/H, ..., alpha/H)
public:
  //---- Hyperparameters ----//
  // number of classes
  double H0;
  // vector of dirichlet parameters
  arma::vec alpha0;
  // Number of observations
  double n;
  // Allocation vector
  arma::uvec C;

  std::random_device rd;
public:
  //=========== CONSTRUCTORS =============//
  // Defaults constructor
  priorDirichlet(double H, double n_obs);

  // Same parameter for each class
  priorDirichlet(double alpha, double H, double n_obs);
  // pass vector of parameters
  priorDirichlet(arma::vec alpha, double n_obs);
  //=========== METHODS =============//
  // Update alpha vector
  arma::vec updateWeigths(arma::vec n_classes);
  // Sample categorical vector
  arma::uvec sampleClass(arma::vec nu);
  // Log conditional distrbution
  arma::vec log_conditional(arma::vec n_class);
};

priorDirichlet::priorDirichlet(double H, double n_obs)
{
  H0 = H;
  n = n_obs; C.set_size(n);
  alpha0.set_size(H); alpha0.fill(1/H);
}

//=========== CONSTRUCTORS =============//
priorDirichlet::priorDirichlet(double alpha, double H, double n_obs)
{
  alpha0.set_size(H); alpha0.fill(alpha/H);
  n = n_obs; H0 = H;
  C.set_size(n);
}

priorDirichlet::priorDirichlet(arma::vec alpha, double n_obs)
{
  H0 = alpha.n_elem;
  alpha0 = alpha;
  n = n_obs;
  C.set_size(n);
}

//=========== METHODS =============//
// Update dirichlet weights - for each i - alpha_i + n_i
arma::vec priorDirichlet::updateWeigths(arma::vec n_class)
{
  NumericMatrix nu(H0);
  arma::vec alpha_star = alpha0 + n_class;
  nu = cpp_rdirichlet(1, as<NumericMatrix>(wrap(alpha_star.t())));
  return as<arma::vec>(wrap(nu));
}

//-------------------------------//
// Sample from discrete distribution with probability vector \nu
// Output is vector of length n

arma::uvec priorDirichlet::sampleClass(arma::vec nu)
{
  NumericVector labels = as<NumericVector>(wrap(seq_len(H0) - 1)); // class from 0 to H-1

  for (int i = 0; i < n; i++)
  {
    C(i) = as<int>(wrap(Rcpp::RcppArmadillo::sample(labels, 1 , TRUE, as<NumericVector>(wrap(nu.t())))));
  }

  // output vector of integer values
  return C;
}
//-------------------------------//
// Conditional distribution in log scale for obs n + 1

arma::vec priorDirichlet::log_conditional( arma::vec n_class)
{
  arma::vec log_conditional(H0);
  for(int h = 0; h < H0; h++)
  {
    log_conditional(h) = log(n_class(h) + alpha0(h)) - log(sum(alpha0) + n - 1);
  }
  return log_conditional;
}

#endif
