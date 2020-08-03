#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// using namespace arma;

double entropy(arma::uvec cl) {

  // This function compute the entropy of a clustering

  // Number of observations
  double n = cl.n_elem;
  // Clustering entropy
  double H_c = 0;
  //--------------------------------------//
  // Tables useful for the computation
  arma::vec counts_set(n, arma::fill::zeros);

  for(int i = 0; i < n ; i++)
  {
    // Compute block sizes for the first clustering
    counts_set(cl[i]) += 1.0;
  }
  // Compute proportion for non-zero elements
  counts_set(find(counts_set > 0))= counts_set(find(counts_set > 0))/n;

  // Compute entropy of clustering
  H_c = -sum(nonzeros(counts_set)%log2(nonzeros(counts_set)));

  return  H_c;
}

// [[Rcpp::export]]
double vi_distC(arma::uvec cl1, arma::uvec cl2) {

  // This function compute the variation of information
  // between two clusterings

  // Number of observations
  int n = cl2.n_elem;

  // Clustering entropy
  double H_c1 = 0;
  double H_c2 = 0;
  // Clusterings mutual information
  double H_cc = 0;

  //--------------------------------------//
  // Tables useful for the computation
  arma::mat counts_int(n,n, arma::fill::zeros);

  for(int i = 0; i < n ; i++)
  {
    // Compute sizes of intersection
    counts_int(cl1[i], cl2[i]) += 1.0;
    // Rcout << i << " " <<  counts_int(cl1[i], cl2[i]) << "\n";
  }
  // Compute proportions for non-zero cells
  counts_int(find(counts_int > 0)) = counts_int(find(counts_int > 0))/n ;

  // Compute entropy of clusterings
  H_c1 = entropy(cl1);
  H_c2 = entropy(cl2);
  // Compute mutual information
  H_cc = -(sum(nonzeros(counts_int)%log2(nonzeros(counts_int))));

  // VI distance
  return  (2*H_cc - H_c1 - H_c2);

}
