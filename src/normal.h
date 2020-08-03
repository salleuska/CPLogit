// PolyaGamma Logistic regression
#ifndef PGLOGIT_HPP
#define PGLOGIT_HPP

#include <RcppArmadillo.h>
#include "marginal_vi.h"
#include "rndpp_mvnormal.h"
#include "PG.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

class PGLOGIT{
  // Polya gamma logistic regression logit(\pi) = \alpha + x_i^^T \beta
  // @n_iter = number of iterations
  // @Y      = response variable
  // @X = explanatory variables
  // @p      = number of variables

public:
  // ---- DIMENSIONS ---- //
  // Number of observations,, variables
  double n; double p;
  // number of iterations
  int n_iter;
  // Data - regression matrix and response vector
  arma::mat X;
  arma::vec Y;
  // ---- regression paramters ---- //
  // @intercept  = intercept
  // @beta       = regression coefficients
  // @omega      = polya-gamma latent variables
  arma::vec intercept; arma::mat beta;
  arma::vec omega;
  // ---- Hypermarameters ---- //
  // Intercept ~ N(a0, tau0)
  double a0; double tau0;
  // Beta coefficients ~ Np(b0, Q)
  arma::vec b0; arma::mat Q0;
  // ---- utilities ---- //

public:
  //=========== CONSTRUCTOR =============//
  // Create the model and allocate objects related to quntities of interest
  PGLOGIT(arma::vec Y0, arma::mat X0, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0);

  //=========== METHODS =============//
  // Initialize prior values
  void setPrior();
  void updateIntercept(int iter);
  void updatePolyaVars(int iter);
  void updateBeta(int iter);
};


//=========== CONSTRUCTOR =============//
PGLOGIT::PGLOGIT(arma::vec Y0, arma::mat X0, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0)
{
  n_iter = n_iter0;

  X = X0; Y = Y0;
  // Number of observations and variables
  n = X0.n_rows; p = X0.n_cols;
  // Set Hypermarameters
  a0 = a_prior; tau0 = tau_prior;
  b0 = b_prior; Q0 = Q_prior;
  // Set size of parameters of interest
  intercept.zeros(n_iter);
  beta.zeros(p, n_iter);
  // Latent polya-gamma variables
  omega.zeros(n);
}
//=========== METHODS =============//
//---------------------------------//
// Set prior values
//---------------------------------//
void PGLOGIT::setPrior()
{
    // Sample intercept from prior distribution
    intercept(0) = ::Rf_rnorm(a0, sqrt(1/tau0));
    // Sample coefficients from multivariate normal distributions
    beta.col(0) = rndpp_mvnormal(1, b0, Q0).t();
}

//------------------------------------//
// update Polya-Gamma latent variables
//------------------------------------//

void PGLOGIT::updatePolyaVars(int iter)
{
   arma::vec pgshape(n); arma::vec pgscale(n);
   pgshape.ones();
   pgscale = intercept(iter - 1) + X*beta.col(iter -1);
   omega = rpg(pgshape, pgscale);
   omega = rpg(pgshape, pgscale);
}

//---------------------------------//
// update intercept term
//---------------------------------//
void PGLOGIT::updateIntercept(int iter)
{
    // Update precision parameter
   double tau_star= tau0 + sum(omega);
   // Compute k values using polya-gamma augmented vars and update mean parameter
   arma::vec xbeta = (X*beta.col(iter -1));
   double k = sum(Y - 0.5 - omega%xbeta);
   double a_star = (a0*tau0 + k)/tau_star;
   intercept(iter) = ::Rf_rnorm(a_star, sqrt(1/tau_star));
}
//------------------------------------//
// update coefficients for each class H
//------------------------------------//
void PGLOGIT::updateBeta(int iter)
{
  // Useful quantities to compute the updated parameters
  arma::mat quadratic_form(p, p, arma::fill::zeros);
  arma::mat xtk_tmp(p, 1, arma::fill::zeros);
  arma::uvec sel;

  // utilities
  arma::mat half_prod;
  arma::vec k_j;
  // updated parameters
  arma::mat Q_star(p, p, arma::fill::zeros);
  arma::vec b_star(p, arma::fill::zeros);

  // compute quadratic form X^T Omega X
  half_prod = X.each_col()%sqrt(omega);
  quadratic_form = half_prod.t()*half_prod;

  // compute polya-gamma residuals
  k_j = Y - 0.5 - omega*intercept(iter);
  xtk_tmp = X.t()*k_j;

  Q_star = arma::inv(quadratic_form + arma::inv(Q0));
  b_star = Q_star*(xtk_tmp + arma::inv(Q0)*b0);

  // sample value for beta
  beta.col(iter) = rndpp_mvnormal(1, b_star, Q_star).t();

}


#endif
