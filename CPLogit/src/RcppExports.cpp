// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dist_from
List dist_from(arma::uvec clust, bool return_partitions);
RcppExport SEXP _CPLogit_dist_from(SEXP clustSEXP, SEXP return_partitionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type clust(clustSEXP);
    Rcpp::traits::input_parameter< bool >::type return_partitions(return_partitionsSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_from(clust, return_partitions));
    return rcpp_result_gen;
END_RCPP
}
// DPgibbsLinear
List DPgibbsLinear(List X_list0, List Y_list0, arma::uvec Z_prior, double psi_par, double H, double H_upper, double alpha, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0);
RcppExport SEXP _CPLogit_DPgibbsLinear(SEXP X_list0SEXP, SEXP Y_list0SEXP, SEXP Z_priorSEXP, SEXP psi_parSEXP, SEXP HSEXP, SEXP H_upperSEXP, SEXP alphaSEXP, SEXP a_priorSEXP, SEXP tau_priorSEXP, SEXP b_priorSEXP, SEXP Q_priorSEXP, SEXP n_iter0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X_list0(X_list0SEXP);
    Rcpp::traits::input_parameter< List >::type Y_list0(Y_list0SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type Z_prior(Z_priorSEXP);
    Rcpp::traits::input_parameter< double >::type psi_par(psi_parSEXP);
    Rcpp::traits::input_parameter< double >::type H(HSEXP);
    Rcpp::traits::input_parameter< double >::type H_upper(H_upperSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type a_prior(a_priorSEXP);
    Rcpp::traits::input_parameter< double >::type tau_prior(tau_priorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b_prior(b_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_prior(Q_priorSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter0(n_iter0SEXP);
    rcpp_result_gen = Rcpp::wrap(DPgibbsLinear(X_list0, Y_list0, Z_prior, psi_par, H, H_upper, alpha, a_prior, tau_prior, b_prior, Q_prior, n_iter0));
    return rcpp_result_gen;
END_RCPP
}
// DPgibbsLogit
List DPgibbsLogit(List X_list0, List Y_list0, arma::uvec Z_prior, double psi_par, double H, double alpha, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0);
RcppExport SEXP _CPLogit_DPgibbsLogit(SEXP X_list0SEXP, SEXP Y_list0SEXP, SEXP Z_priorSEXP, SEXP psi_parSEXP, SEXP HSEXP, SEXP alphaSEXP, SEXP a_priorSEXP, SEXP tau_priorSEXP, SEXP b_priorSEXP, SEXP Q_priorSEXP, SEXP n_iter0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X_list0(X_list0SEXP);
    Rcpp::traits::input_parameter< List >::type Y_list0(Y_list0SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type Z_prior(Z_priorSEXP);
    Rcpp::traits::input_parameter< double >::type psi_par(psi_parSEXP);
    Rcpp::traits::input_parameter< double >::type H(HSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type a_prior(a_priorSEXP);
    Rcpp::traits::input_parameter< double >::type tau_prior(tau_priorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b_prior(b_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_prior(Q_priorSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter0(n_iter0SEXP);
    rcpp_result_gen = Rcpp::wrap(DPgibbsLogit(X_list0, Y_list0, Z_prior, psi_par, H, alpha, a_prior, tau_prior, b_prior, Q_prior, n_iter0));
    return rcpp_result_gen;
END_RCPP
}
// gibbsPGLogit
List gibbsPGLogit(arma::vec Y0, arma::mat X0, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0);
RcppExport SEXP _CPLogit_gibbsPGLogit(SEXP Y0SEXP, SEXP X0SEXP, SEXP a_priorSEXP, SEXP tau_priorSEXP, SEXP b_priorSEXP, SEXP Q_priorSEXP, SEXP n_iter0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y0(Y0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< double >::type a_prior(a_priorSEXP);
    Rcpp::traits::input_parameter< double >::type tau_prior(tau_priorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b_prior(b_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_prior(Q_priorSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter0(n_iter0SEXP);
    rcpp_result_gen = Rcpp::wrap(gibbsPGLogit(Y0, X0, a_prior, tau_prior, b_prior, Q_prior, n_iter0));
    return rcpp_result_gen;
END_RCPP
}
// gibbsLinearModel
List gibbsLinearModel(arma::vec Y0, arma::mat X0, double a_prior, double tau_prior, arma::vec b_prior, arma::mat Q_prior, int n_iter0);
RcppExport SEXP _CPLogit_gibbsLinearModel(SEXP Y0SEXP, SEXP X0SEXP, SEXP a_priorSEXP, SEXP tau_priorSEXP, SEXP b_priorSEXP, SEXP Q_priorSEXP, SEXP n_iter0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y0(Y0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< double >::type a_prior(a_priorSEXP);
    Rcpp::traits::input_parameter< double >::type tau_prior(tau_priorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b_prior(b_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_prior(Q_priorSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter0(n_iter0SEXP);
    rcpp_result_gen = Rcpp::wrap(gibbsLinearModel(Y0, X0, a_prior, tau_prior, b_prior, Q_prior, n_iter0));
    return rcpp_result_gen;
END_RCPP
}
// rndpp_mvnormal
arma::mat rndpp_mvnormal(int n, const arma::vec& mean, const arma::mat& sigma);
RcppExport SEXP _CPLogit_rndpp_mvnormal(SEXP nSEXP, SEXP meanSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rndpp_mvnormal(n, mean, sigma));
    return rcpp_result_gen;
END_RCPP
}
// vi_distC
double vi_distC(arma::uvec cl1, arma::uvec cl2);
RcppExport SEXP _CPLogit_vi_distC(SEXP cl1SEXP, SEXP cl2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type cl1(cl1SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type cl2(cl2SEXP);
    rcpp_result_gen = Rcpp::wrap(vi_distC(cl1, cl2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CPLogit_dist_from", (DL_FUNC) &_CPLogit_dist_from, 2},
    {"_CPLogit_DPgibbsLinear", (DL_FUNC) &_CPLogit_DPgibbsLinear, 12},
    {"_CPLogit_DPgibbsLogit", (DL_FUNC) &_CPLogit_DPgibbsLogit, 11},
    {"_CPLogit_gibbsPGLogit", (DL_FUNC) &_CPLogit_gibbsPGLogit, 7},
    {"_CPLogit_gibbsLinearModel", (DL_FUNC) &_CPLogit_gibbsLinearModel, 7},
    {"_CPLogit_rndpp_mvnormal", (DL_FUNC) &_CPLogit_rndpp_mvnormal, 3},
    {"_CPLogit_vi_distC", (DL_FUNC) &_CPLogit_vi_distC, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_CPLogit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
