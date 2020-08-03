// Function to test GLogit class

// uff#include "gibbs.h"
#include "GLogit.h"
#include "GLinearModel.h"
#include "PGLogistic.h"
#include "DPGLogit.h"
#include "DPGLinear.h"
#include "linearModel.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]
List DPgibbsLinear(List X_list0, List Y_list0, arma::uvec Z_prior, double psi_par,  double H, double H_upper, double alpha,double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0){
  
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


// [[Rcpp::export]]
List gibbs(List X_list0, List Y_list0, arma::uvec Z_prior, double psi_par,  double H, double alpha,double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior,  int n_iter0){

  GLOGIT my_mod(X_list0, Y_list0, Z_prior, psi_par,  H, alpha, a_prior, tau_prior, b_prior, Q_prior, n_iter0);

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

// [[Rcpp::export]]
List DPgibbs(List X_list0, List Y_list0, arma::uvec Z_prior, double psi_par,  double H, double H_upper, double alpha,double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0){

  DPGLOGIT my_mod(X_list0, Y_list0, Z_prior, psi_par, H, H_upper, alpha, a_prior, tau_prior, b_prior, Q_prior, n_iter0);

  Rcout << "set prior function\n";
  my_mod.setPrior();

  Rcout << "Start gibbs\n";

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


// [[Rcpp::export]]
List gibbsPGLogit(arma::vec Y0, arma::mat X0, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0)
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


// [[Rcpp::export]]
List gibbsLinearModel(arma::vec Y0, arma::mat X0, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0)
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

// [[Rcpp::export]]
List gibbsGLinearModel(List X_list0, List Y_list0, arma::uvec Z_prior, double psi_par,  double H, double alpha,double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior,  int n_iter0){

  GLINEAR my_mod(X_list0, Y_list0, Z_prior, psi_par,  H, alpha, a_prior, tau_prior, b_prior, Q_prior, n_iter0);

  Rcout << "setPriorfunction" << "\n";
  my_mod.setPrior();

  for(int t = 1; t < n_iter0; ++t)
  {
    // my_mod.updatePolyaVars(t);
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

