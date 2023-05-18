// Function to test GLogit class

// uff#include "gibbs.h"
#include "GLogit.h"
#include "GLinearModel.h"
#include "PGLogistic.h"
#include "DPGLogit.h"
#include "DPGLogitMultiplePrior.h"
#include "DPGLinear.h"
#include "DPGLinearNoInt.h"
#include "linearModel.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Grouped linear regression using CP Process prior with DP base EPPF.
//'
//' This function does a gibbs sampling for a linear regression model using a CP prior for 
//' clustering regression coefficients. 
//'
//' @param X_list0 list of n_i x p matrices of covariates (including intercept)
//' @param Y_list0 list of n_i x 1 vectors of response variable
//' @param Z_prior prior information on clustering
//' @param psi_par value for the penalizations parameter
//' @param H number or clusters
//' @param H_upper upper limit for number of clusters
//' @param alpha value for the Dirichlet Process concentration parameter
//' @param a_prior mean of the normal prior on intercept
//' @param tau_prior inverse of variance of the normal prior on intercept
//' @param b_prior mean vector of the normal multivariate prior for the coefficients
//' @param Q_prior covariance matrix of the normal multivariate prior for the coefficients
//' @param n_iter0 total number of iterations
//' @export
// [[Rcpp::export]]
List gibbsLinearCP(List X_list0, List Y_list0, arma::uvec Z_prior, double psi_par,  double H, double H_upper, double alpha,double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0){
  
  DPGLINEAR my_mod(X_list0, Y_list0, Z_prior, psi_par, H, H_upper, alpha, a_prior, tau_prior, b_prior, Q_prior, n_iter0);
  
  Rcout << "set prior function\n";
  my_mod.setPrior();
  
  Rcout << "Start gibbs\n";
  
  for(int t = 1; t < n_iter0; ++t)
  {
    // Rcout << "Update intercept\n";
    my_mod.updateIntercept(t);
    // Rcout << "update beta\n";
    my_mod.updateBeta(t);
    // Rcout << "update class\n";
    my_mod.updateClass(t);
    
    if (t % 100 == 0) Rcout << t << "\n";
    
  }
  //--------- End Gibbs --------//
  // return estimated quantities;
  return List::create(_["clustering"] = my_mod.clustering,
                      _["intercept"] = my_mod.intercept,_["beta"] = my_mod.beta);
  
}


//' Grouped logistic regression using Centered partition process.
//'
//' This function does a gibbs sampling for a logistic regression model using a DP prior for 
//' clustering regression coefficients. Sampling for the logistic regressions uses Polya-Gamma 
//' data augmentation.
//'
//' @param X_list0 list of n_i x p matrices of covariates (including intercept)
//' @param Y_list0 list of n_i x 1 vectors of response variable
//' @param Z_prior initial values for clustering allocation
//' @param psi_par value for the penalizations parameter
//' @param H number or clusters
//' @param H_upper upper limit for number of clusters
//' @param alpha value for the Dirichlet Process concentration parameter
//' @param a_prior mean of the normal prior on intercept
//' @param tau_prior inverse of variance of the normal prior on intercept
//' @param b_prior mean vector of the normal multivariate prior for the coefficients
//' @param Q_prior covariance matrix of the normal multivariate prior for the coefficients
//' @param n_iter0 total number of iterations
//' @export
// [[Rcpp::export]]
List gibbsLogitCP(List X_list0, List Y_list0, arma::uvec Z_prior, double psi_par,  double H, double H_upper,
  double alpha,double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior,  int n_iter0){

  DPGLOGIT my_mod(X_list0, Y_list0, Z_prior, psi_par,  H, H_upper, alpha, a_prior, tau_prior, b_prior, Q_prior, n_iter0);

  Rcout << "setPriorfunction" << "\n";
  my_mod.setPrior();

  for(int t = 1; t < n_iter0; ++t)
  {
    my_mod.updatePolyaVars(t);
    my_mod.updateIntercept(t);
    my_mod.updateBeta(t);
    my_mod.updateClass(t);

    if (t % 100 == 0) Rcout << t << "\n";
  }
  //--------- End Gibbs --------//
  // return estimated quantities;
  return List::create(_["clustering"] = my_mod.clustering,
                    _["intercept"] = my_mod.intercept,_["beta"] = my_mod.beta);

}

//' Grouped logistic regression using Centered Partition process - multiple prior partitions
//'
//' This function does a gibbs sampling for a logistic regression model using a DP prior for 
//' clustering regression coefficients. Sampling for the logistic regressions uses Polya-Gamma 
//' data augmentation.
//'
//' @param X_list0 list of n_i x p matrices of covariates (including intercept)
//' @param Y_list0 list of n_i x 1 vectors of response variable
//' @param Z_prior matrix of possible partitions for clustering allocation
//' @param psi_par value for the penalizations parameter
//' @param H number or clusters
//' @param H_upper upper limit for number of clusters
//' @param alpha value for the Dirichlet Process concentration parameter
//' @param a_prior mean of the normal prior on intercept
//' @param tau_prior inverse of variance of the normal prior on intercept
//' @param b_prior mean vector of the normal multivariate prior for the coefficients
//' @param Q_prior covariance matrix of the normal multivariate prior for the coefficients
//' @param n_iter0 total number of iterations
//' @export
// [[Rcpp::export]]
List gibbsLogitCPMultiple(List X_list0, List Y_list0, arma::umat Z_prior, double psi_par,  double H, double H_upper,
  double alpha,double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior,  int n_iter0){

  DPGLOGITMULTIPLE my_mod(X_list0, Y_list0, Z_prior, psi_par,  H, H_upper, alpha, a_prior, tau_prior, b_prior, Q_prior, n_iter0);

  Rcout << "setPriorfunction" << "\n";
  my_mod.setPrior();
  

  for(int t = 1; t < n_iter0; ++t)
  {
    my_mod.updatePolyaVars(t);
    my_mod.updateIntercept(t);
    my_mod.updateBeta(t);
    my_mod.updateClass(t);

    if (t % 100 == 0) Rcout << t << "\n";
  }
  //--------- End Gibbs --------//
  // return estimated quantities;
  return List::create(_["clustering"] = my_mod.clustering,
                    _["intercept"] = my_mod.intercept,_["beta"] = my_mod.beta);

}

//' Logistic regression using Polya-Gamma data augmentation.
//'
//' This function does a gibbs sampling for a logistic regression using Polya-Gamma 
//' data augmentation.
//'
//' @param Y n x 1 vector of response variable 
//' @param X n x p matrix of covariates (including intercept)
//' @param tau_prior inverse of variance of the normal prior on intercept
//' @param a_prior mean of the normal prior on intercept
//' @param b_prior mean vector of the normal multivariate prior for the coefficients
//' @param Q_prior covariance matrix of the normal multivariate prior for the coefficients
//' @param n_iter0 total number of iterations
//' @export
// [[Rcpp::export]]
List gibbsLogit(arma::vec Y0, arma::mat X0, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0)
{

  PGLOGIT my_mod( Y0, X0, a_prior, tau_prior,  b_prior,  Q_prior, n_iter0);

  Rcout << "setPriorfunction" << "\n";
  my_mod.setPrior();

  for(int t = 1; t < n_iter0; ++t)
  {
    my_mod.updatePolyaVars(t);
    my_mod.updateIntercept(t);

    my_mod.updateBeta(t);

    if (t % 100 == 0) Rcout << t << "\n";
  }
  //--------- End Gibbs --------//
  // return estimated quantities;
  return List::create(_["intercept"] = my_mod.intercept,_["beta"] = my_mod.beta);

}


//' Bayesian linear regression.
//'
//' This function does a gibbs sampling for a linear regression model.
//'
//' @param Y0 n x 1 vector of response variable 
//' @param X0 n x p matrix of covariates (including intercept)
//' @param tau_prior inverse of variance of the normal prior on intercept
//' @param a_prior mean of the normal prior on intercept
//' @param b_prior mean vector of the normal multivariate prior for the coefficients
//' @param Q_prior covariance matrix of the normal multivariate prior for the coefficients
//' @param n_iter0 total number of iterations
//' @export
// [[Rcpp::export]]
List gibbsLinear(arma::vec Y0, arma::mat X0, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0)
{

  LINEARMODEL my_mod( Y0, X0, a_prior, tau_prior,  b_prior,  Q_prior, n_iter0);

  Rcout << "setPriorfunction" << "\n";
  my_mod.setPrior();

  for(int t = 1; t < n_iter0; ++t)
  {
    //    my_mod.updatePolyaVars(t);
    my_mod.updateIntercept(t);
    my_mod.updateBeta(t);
    //    my_mod.updateClass(t);

    if (t % 100 == 0) Rcout << t << "\n";
  }
  //--------- End Gibbs --------//
  // return estimated quantities;
  return List::create(_["intercept"] = my_mod.intercept,_["beta"] = my_mod.beta);

}

//' Grouped Bayesian linear regression - only regressors.
//'
//' This function does a gibbs sampling for a linear regression model with only regressors coef
//' note: this is a trick to clustering means of a normal distributed data
//'
//' @param X_list0 list of n_i x p matrices of covariates (including intercept)
//' @param Y_list0 list of n_i x 1 vectors of response variable
//' @param Z_prior prior information on clustering
//' @param psi_par value for the penalizations parameter
//' @param H number or clusters
//' @param H_upper upper limit for number of clusters
//' @param alpha value for the Dirichlet Process concentration parameter
//' @param b_prior mean vector of the normal multivariate prior for the coefficients
//' @param Q_prior covariance matrix of the normal multivariate prior for the coefficients
//' @param n_iter0 total number of iterations
//' @export
// [[Rcpp::export]]
List gibbsDensityCP(List X_list0, List Y_list0, arma::uvec Z_prior, double psi_par,  double H, double H_upper, double alpha, arma::vec b_prior, arma::mat Q_prior, int n_iter0){
  
  DPGLINEARNOINT my_mod(X_list0, Y_list0, Z_prior, psi_par, H, H_upper, alpha, b_prior, Q_prior, n_iter0);
  
  Rcout << "set prior function\n";
  my_mod.setPrior();
  
  Rcout << "Start gibbs\n";
  
  for(int t = 1; t < n_iter0; ++t)
  {
    // Rcout << "Update intercept\n";
    // my_mod.updateIntercept(t);
    // Rcout << "update beta\n";
    my_mod.updateBeta(t);
    // Rcout << "update class\n";
    my_mod.updateClass(t);
    
    if (t % 100 == 0) Rcout << t << "\n";
    
  }
  //--------- End Gibbs --------//
  // return estimated quantities;
  return List::create(_["clustering"] = my_mod.clustering,
                      _["beta"] = my_mod.beta);
  
}

