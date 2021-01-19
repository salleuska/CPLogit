// Grouped Linear regression for multiple dataset
#ifndef DPGLINEAR_H
#define DPGLINEAR_H

#include <RcppArmadillo.h>
#include "priorDirichlet.h"
#include "marginal_vi.h"
#include "rndpp_mvnormal.h"
// #include "PG.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

class DPGLINEAR{
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
  // Upper limit for number of classes
  double H_up;
  // Vector of class labels
  IntegerVector class_labels;
  // List of datasets - regression matrix and response vector
  arma::field<arma::mat> X_list;
  arma::field<arma::vec> Y_list;
  // base clustering
  arma::uvec Z0; double psi; int H_Z0;
  // ---- PRIOR DISTRIBUTION ---- //
  // @prior      = Prior distribution for mixture classes
  priorDP prior;
  // ---- regression parameters ---- //
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
  DPGLINEAR(List X_list0, List Y_list0, arma::uvec Z_prior, double psi_par,  double H, double H_upper,  double alpha, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0);
  
  //=========== METHODS =============//
  // Initialize prior values
  void setPrior();
  void updateIntercept(int iter);
  void updateBeta(int iter);
  void updateClass(int iter);
};

//=========== CONSTRUCTOR =============//
DPGLINEAR::DPGLINEAR(List X_list0, List Y_list0, arma::uvec Z_prior,double psi_par, double H, double H_upper, double alpha, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0) : prior(H, H_upper, alpha, X_list0.size())
{
  // Number of datastes
  n = X_list0.size();
  // Base partition
  Z0 = Z_prior; psi = psi_par;
  // HERE _ ADD SIZE OF clustering base
  arma::uvec Z0_unique = unique(Z0);
  H_Z0 = Z0_unique.n_elem;
  
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
  H0 = H; H_up = H_upper; n_iter = n_iter0;
  // Set Hypermarameters
  a0 = a_prior; tau0 = tau_prior;
  b0 = b_prior; Q0 = Q_prior;
  // Set size of parameters of interest
  intercept.zeros(n, n_iter);
  beta.zeros(H_up, p, n_iter);
  // Cluster allocation and class labels
  clustering.set_size(n, n_iter0);
  class_labels = as<IntegerVector>(wrap(seq_len(H0) - 1));
  // Utils
  class_probs.zeros(n, H_up, n_iter0);
  // class_probs_tmp.zeros(n, n_iter0);
  class_size.zeros(H_up);
}
//=========== METHODS =============//
//---------------------------------//
// Set prior values
//---------------------------------//

void DPGLINEAR::setPrior()
{
  arma::uvec temp_clust(n);
  // Random allocation of observations to classes - label 0 to H0-1
  for (int i = 0; i < n; ++i)
  {
    // temporary labels
    temp_clust(i)= as<int>(wrap(Rcpp::sample(class_labels, 1, FALSE)));
    // sample intercep for each dataset
    intercept(i, 0) = ::Rf_rnorm(a0, sqrt(1/tau0));
  }
  
  // rename labels to have no gaps
  arma::uvec unique_lab = unique(temp_clust);
  if(unique_lab.n_elem < H0)
  {
    clustering.col(0) = relabel(temp_clust);
    H0 = unique_lab.n_elem;
  } else {
    clustering.col(0) = temp_clust;
  }
  // Sample coefficients from multivariate normal distributions
  for(int h = 0; h < H0; h++)
  {
    beta.slice(0).row(h) = rndpp_mvnormal(1, b0, Q0);
  }
}

//---------------------------------//
// update intercept term
//---------------------------------//
void DPGLINEAR::updateIntercept(int iter)
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
void DPGLINEAR::updateBeta(int iter){
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
        xtk_tmp = X_list(sel(j)).t()*Y_list(sel(j));
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
void DPGLINEAR::updateClass(int iter)
{
  // Save current clustering configuration - to update sequentially
  arma::uvec clustering_to_update = clustering.col(iter-1);
  // utilities - clustering part
  arma::uvec unique_vals; int H_minus_i;
  arma::uvec clust_minus_i(n, arma::fill::zeros);
  arma::vec class_size_minus_i(n, arma::fill::zeros);
  NumericVector class_lab_minus_i;
  
  // utilities - likelihood
  arma::vec reg; double log_lik;
  arma::mat beta_new(1,p, arma::fill::zeros);
  // allocation probability vector for obs i
  arma::vec log_p_tmp; arma::vec p_class;
  
  // penalization
  arma::vec dist;
  
  
  for(int i = 0; i < n; ++i)
  {
    // copy current clustering configurations
    clust_minus_i = clustering_to_update;
    // remove i-th observation and compute the number of unique values
    clust_minus_i.shed_row(i);
    unique_vals = unique(clust_minus_i);
    H_minus_i = unique_vals.n_elem;
    
    if(H0 - H_minus_i == 0){
      
      // Observation is NOT a singleton
      log_p_tmp.set_size(H_minus_i + 1);
      dist.set_size(H_minus_i + 1);
      
      
      // Obtain vector of class sizes from the available one
      // by decreasing by one the corresponding class count
      class_size(clustering_to_update(i)) -= 1;
      class_size_minus_i = class_size.head(H_minus_i);
      // Allocation probability to an observed class
      log_p_tmp.head(H_minus_i) = prior.log_conditional(class_size_minus_i);
      
      //-------------------------------------------------
      // Compute lilkelihood for observed classes
      for(int h = 0; h < H_minus_i; ++h)
      {
        reg = (intercept(i, iter) + X_list(i)*beta.slice(iter).row(h).t());
        log_lik = -0.5*sum(arma::square(Y_list(i) - reg));
        log_p_tmp(h) += log_lik;
      }
      //---------------
      // Allocation probability to a new class - Using alg 8 from Neal(2000) m = 1 (lazy)
      // sample a new value for beta for the class H_minus_i + 1
      beta_new = rndpp_mvnormal(1, b0, Q0);
      
      reg = (intercept(i, iter) + X_list(i)*beta_new.t());
      log_p_tmp.tail(1) = prior.log_new_class -0.5*sum(arma::square(Y_list(i) - reg));
      //-------------------------------------------------
      // Add penalization
      if(psi > 0){
        dist = marginal_vi(clustering_to_update, Z0,  H_minus_i + 1,i);
        
        log_p_tmp = log_p_tmp -  psi*dist;
      }
      
      
      //-------------------------------------------------
      //Normalize probabilities
      p_class = normalize_probs(log_p_tmp);
      //-------------------------------------------------
      // Sample class for observation i
      class_lab_minus_i = as<NumericVector>(wrap(seq_len(H_minus_i+1) -1));
      clustering_to_update(i)= as<int>(wrap(Rcpp::sample(class_lab_minus_i, 1, FALSE, as<NumericVector>(wrap(p_class.t())))));
      // If a new class is sampled assign beta_new kernel
      if(clustering_to_update(i) == H_minus_i) {beta.slice(iter).row(H_minus_i) = beta_new;}
      
      
      // Update number of clusters
      unique_vals = unique(clustering_to_update);
      H0 = unique_vals.n_elem;
      // Update class_sizes
      class_size(clustering_to_update(i)) += 1;
      
    } else {
      
      
      // Observation is a singleton
      dist.set_size(H_minus_i + 1);
      log_p_tmp.set_size(H_minus_i + 1);
      
      ///================ REORDERING PART ================== ///
      // Need to fix parameters associated to the cluster of the ith observation
      arma::mat beta_tmp = beta.slice(iter);
      
      beta_tmp.row(clustering_to_update(i)).zeros();
      
      // Reorder parameters such that there are no gaps
      arma::uvec index(H_minus_i + 1);
      // Create ordering vector
      index.head(H_minus_i) = unique_vals;
      index.tail(1) = clustering_to_update(i);
      // Actually reorder vector
      beta_tmp.rows(0, H_minus_i) = beta_tmp.rows(index);
      
      // Substitute in the object
      beta.slice(iter) = beta_tmp;
      //.....................................................//
      // Obtain vector of class sizes from the available one
      // by decreasing by one the corresponding class count
      // Update class sizes removing the i observation redordering
      class_size(clustering_to_update(i)) -= 1;
      class_size.head(H_minus_i +1) = class_size.elem(index);
      
      
      //-------------------------------------------------
      // compute penalization
      // Add penalization
      if(psi > 0){
        dist = marginal_vi(clustering_to_update, Z0, H_minus_i + 1,i);
      }
      
      // Relabel current clustering vector without i obs
      clustering_to_update.shed_row(i);
      clustering_to_update = relabel(clust_minus_i);
      ///================ ============== ================== ///
      // Vector of log probabilties
      log_p_tmp.head(H_minus_i) = prior.log_conditional(class_size.head(H_minus_i));
      
      //-------------------------------------------------
      // Compute lilkelihood for observed classes
      
      for(int h = 0; h < H_minus_i; ++h)
      {
        reg = (intercept(i, iter) + X_list(i)*beta.slice(iter).row(h).t());
        log_lik = -0.5*sum(arma::square(Y_list(i) - reg));
        log_p_tmp(h) += log_lik;
      }
      
      
      //---------------
      // Allocation probability to a new class - Using alg 8 from Neal(2000) m = 1 (lazy)
      // sample a new value for beta for the class H_minus_i + 1
      beta_new = rndpp_mvnormal(1, b0, Q0);
      reg = (intercept(i, iter) + X_list(i)*beta_new.t());
      log_p_tmp.tail(1) = prior.log_new_class -0.5*sum(arma::square(Y_list(i) - reg));
      
      // Check for penalization
      if(psi > 0){ log_p_tmp = log_p_tmp - psi*dist;}
      
      //-------------------------------------------------
      //Normalize probabilities
      p_class = normalize_probs(log_p_tmp);
      
      
      //-------------------------------------------------
      // Sample class for observation i
      class_lab_minus_i = as<NumericVector>(wrap(seq_len(H_minus_i+1) -1));
      int class_for_i = as<int>(wrap(Rcpp::sample(class_lab_minus_i, 1, FALSE, as<NumericVector>(wrap(p_class.t())))));
      
      // If a new class is sampled assign beta_new kernel
      if(class_for_i == H_minus_i) {beta.slice(iter).row(H_minus_i) = beta_new;}
      
      // Update clustering
      clustering_to_update.insert_rows(i, 1);
      clustering_to_update(i) = class_for_i;
      unique_vals = unique(clustering_to_update);
      
      // Update number of clusters
      H0 = unique_vals.n_elem;
      class_size(clustering_to_update(i)) += 1;
    }
    // ---- end ifelse ----- //
    // reset utils vectors
    class_size_minus_i.reset(); H_minus_i = 0;
    clust_minus_i.reset();
    log_p_tmp.reset(); p_class.reset();
    unique_vals.reset(); clust_minus_i.reset();
    dist.reset();
    // class_lab_minus_i.assign(class_lab_minus_i.size(), 0);
    
    
    
    // ---- end for ----- //
  }
  
  
  clustering.col(iter) = clustering_to_update;
  // ---- end function ----- //
}

#endif
