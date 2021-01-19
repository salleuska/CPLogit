// Grouped Linear regression for multiple dataset
#ifndef GLINEAR_H
#define GLINEAR_H

#include <RcppArmadillo.h>
#include "priorDirichlet.h"
#include "marginal_vi.h"
#include "rndpp_mvnormal.h"
// #include "PG.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

class GLINEAR{
  // Grouped logistic regression using exchangeable prior
  // @n_iter = number of iterations
  // @Y_list = list of response variables [j][n_j x 1]
  // @X_list = list of explanatory variables [j][n_j x p]
  // @p      = number of variables
  // @H      = upper bound for number of components

public:
  // ---- DIMENSIONS ---- //
  // Number of datasets, variables
  double n; double p;
  // number of iterations
  int n_iter;
  // Number of elements for each datasets
  arma::uvec n_data;
  // Number of classes
  double H0;
  // Vector of class labels
  IntegerVector class_labels;
  // List of datasets - regression matrix and response vector
  arma::field<arma::mat> X_list;
  arma::field<arma::vec> Y_list;
  // base clustering
  arma::uvec Z0; double psi;
  // ---- PRIOR DISTRIBUTION ---- //
  // @prior     = Prior distribution for mixture classes
    priorDirichlet prior;
  // ---- regression paramters ---- //
  // @intercept  = dataset-specific intercept
  // @beta       = grouped regression coefficients
  // @omega      = polya-gamma latent variables
  // @clustering = clustering vector
  arma::mat intercept; arma::cube beta;
  arma::field<arma::vec> omega;
  arma::umat clustering;
  // ---- Hypermarameters ---- //
  // Intercept ~ N(a0, tau0)
  double a0; double tau0;
  // Beta coefficients ~ Np(b0, Q)
  arma::vec b0; arma::mat Q0;
  // ---- utilities ---- //
  arma::cube class_probs;
  arma::vec class_size;

public:
  //=========== CONSTRUCTOR =============//
  // Create the model and allocate objects related to quntities of interest
  GLINEAR(List X_list0, List Y_list0, arma::uvec Z_prior, double psi_par,  double H, double alpha, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0);

  //=========== METHODS =============//
  // Initialize prior values
  void setPrior();
  void updatePolyaVars(int iter);
  void updateIntercept(int iter);
  void updateBeta(int iter);
  void updateClass(int iter);
};


//=========== CONSTRUCTOR =============//
GLINEAR::GLINEAR(List X_list0, List Y_list0, arma::uvec Z_prior,double psi_par, double H, double alpha, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0) : prior(alpha, H, X_list0.size())
{
  // Number of datastes
  n = X_list0.size();
  // Base partition
  Z0 = Z_prior; psi = psi_par;
  // Save field of data list and populate vector of sizes
  n_data.set_size(n); n_data.fill(0);
  X_list.set_size(n); Y_list.set_size(n);
  for(int k = 0; k < n; ++k)
  {
    // data matrix(n_j x p)
    arma::mat tmp = X_list0[k];
    X_list(k) = tmp;
    arma::vec tmp2 = Y_list0[k];
    Y_list(k) = tmp2;
    n_data(k) = tmp.n_rows;
  }
  p = X_list(0).n_cols;

  // Set upper bound for number of classes and interations number
  H0 = H; n_iter = n_iter0;
  // Set Hypermarameters
  a0 = a_prior; tau0 = tau_prior;
  b0 = b_prior; Q0 = Q_prior;
  // Set size of parameters of interest
  intercept.zeros(n, n_iter);
  beta.zeros(H, p, n_iter);
  // Latent polya-gamma variables
  // omega.set_size(n);
  // for(int j = 0; j < n; ++j)
  // {
  //   arma::vec tmp(n_data(j));
  //   omega(j) = tmp;
  // }
  // Cluster allocation and class labels
  clustering.set_size(n, n_iter0);
  class_labels = as<IntegerVector>(wrap(seq_len(H0) - 1));
  // Utils
  class_probs.zeros(n, H0, n_iter0);
  // class_probs_tmp.zeros(n, n_iter0);
  class_size.zeros(H0);
}
//=========== METHODS =============//
//---------------------------------//
// Set prior values
//---------------------------------//
void GLINEAR::setPrior()
{
  // Random allocation of observations to classes
  for (int i = 0; i < n; ++i)
  {

    clustering(i,0)= as<int>(wrap(Rcpp::RcppArmadillo::sample(class_labels, 1, FALSE)));
    // Sample intercept from prior normal distribution
    intercept(i, 0) = ::Rf_rnorm(a0, sqrt(1/tau0));
  }
  // Sample coefficients from multivariate normal distributions
  for(int h = 0; h < H0; h++)
  {
    beta.slice(0).row(h) = rndpp_mvnormal(1, b0, Q0);
  }
}

//---------------------------------//
// update Polya-Gamma latent variables
//---------------------------------//

// void GLINEAR::updatePolyaVars(int iter)
// {
//  for(int j = 0; j < n; ++j)
//  {
//    arma::uword index_beta = clustering(j, iter -1);
//    arma::vec pgshape(n_data(j)); arma::vec pgscale(n_data(j));
//    pgshape.ones();
//    pgscale = intercept(j, iter - 1) + X_list(j)*beta.slice(iter -1).row(index_beta).t();
//    omega(j) = rpg(pgshape, pgscale);
//  }
// }
//---------------------------------//
// update intercept term
//---------------------------------//
void GLINEAR::updateIntercept(int iter)
{
 for(int j = 0; j < n; ++j)
 {
   // Update precision parameter
   double tau_star= tau0 + n;
   // update mean parameter
   double den = (1/tau0)/(1/tau0 + 1/n);
   double y_mean = sum(Y_list(j))/n;
   double a_star = y_mean*(1/tau0)/den + a0/den;
   intercept(iter) = ::Rf_rnorm(a_star, sqrt(1/tau_star));
   
 }
}
//------------------------------------//
// update coefficients for each class H
//------------------------------------//
void GLINEAR::updateBeta(int iter)
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

  for(int h = 0; h < H0; ++h)
  {
    arma::uvec sel = find(clustering.col(iter-1) == h);
    class_size(h) = sel.n_elem;

    if(class_size(h) > 0){
      for(int j = 0; j < class_size(h); ++j)
      {
        // compute quadratic form X^T  X
        quadratic_form = X_list(sel(j)).t()*X_list(sel(j));
        xtk_tmp = X_list(sel(j))*Y_list(sel(j));
      }

      Q_star = arma::inv(quadratic_form + arma::inv(Q0));
      b_star = Q_star*xtk_tmp;

      // sample value for beta
      beta.slice(iter).row(h) = rndpp_mvnormal(1, b_star, Q_star);

    }
    else
    {
      // sample from the prior distribution
      beta.slice(iter).row(h) = rndpp_mvnormal(1, b0, Q0);
    }

    // reset tmp quantities
    quadratic_form.zeros(); xtk_tmp.zeros();
  }

}

//-----------------------------------------------------------/
// Update class allocation using prior info about clustering
//-----------------------------------------------------------/
void GLINEAR::updateClass(int iter)
{
  // Save current clustering configuration - to update sequentially
  arma::uvec clustering_to_update = clustering.col(iter-1);
  // utilities
  arma::uvec sel; arma::vec reg; double log_lik;
  arma::vec class_size_minus_i(H0, arma::fill::zeros);
  arma::vec prob_tmp(H0, arma::fill::zeros);
  arma::vec dist(H0, arma::fill::zeros);
  arma::vec prob_h(H0, arma::fill::zeros);
  
  for(int j = 0; j < n; ++j)
  {
    arma::vec var(n_data(j), arma::fill::ones);
     for(int h = 0; h < H0; ++h)
     {
       // cluster size removing ith observation
       arma::uvec sel = find(clustering_to_update == h);
       class_size_minus_i(h) = sel.n_elem;
     }
     class_size_minus_i(clustering_to_update(j)) =  class_size_minus_i(clustering_to_update(j)) - 1;
     // Compute class probabilities - log
     prob_tmp = log(class_size_minus_i + prior.alpha0) - log(sum(prior.alpha0) + n - 1);

     for(int h = 0; h < H0; ++h)
     {
       reg = (intercept(j, iter) + X_list(j)*beta.slice(iter).row(h).t());
       log_lik = -0.5*sum(square(Y_list(j) - reg));
       prob_tmp(h) += log_lik;
     }

     dist = marginal_vi( clustering_to_update, Z0, H0, j);
     prob_tmp  = prob_tmp - psi*dist;

     // Compute class probabilities - using the trick exp{-log(1 + sum(...))}
     class_probs.slice(iter).row(j) =  normalize_probs(prob_tmp).t();
     clustering_to_update(j) = as<int>(wrap(Rcpp::RcppArmadillo::sample(class_labels, 1, TRUE, as<NumericVector>(wrap(class_probs.slice(iter).row(j).t())))));
  }

  clustering.col(iter) = clustering_to_update;
}


#endif
