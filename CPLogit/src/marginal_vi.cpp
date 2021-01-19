#include "marginal_vi.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec marginal_vi(arma::uvec clustering, arma::uvec C_0, int H, int obs) {

  // This function compute the "marginal" variation of information
  // between the base clustering C_0 and all possible clusterings that
  // can be created by moving the ith observation

  //@clustering = current clustering
  //@C_0 = base clustering (known)
  //@H = possible classes
  //@obs = obervation


  // coping elements (because you originally were recoding the labels from R)
  arma::uvec clust_minus = clustering;
  arma::uvec C_0_minus = C_0;
  obs = obs;

  // number elements cluster c0
  arma::uvec unique_v = unique(C_0);
  int H_C0 = unique_v.n_elem;

  // Number of observations
  double n = C_0.size();

  // Clustering entropy
  double H_c = 0;
  // Clusterings mutual information
  double H_cc = 0;

  // Vector of distances
  arma::vec dist(H, arma::fill::zeros);

  //--------------------------------------//
  // Tables useful for the computation
  arma::mat counts_int(H,H_C0, arma::fill::zeros);
  arma::vec counts_set(H, arma::fill::zeros);

  for(int i = 0; i < n; i++)
  {
    if(i == obs)
    {
      // Skip subject
      continue;
    }
    else
    {
      // Compute block sizes for the current clustering minus the observation obs
      counts_set(clust_minus[i]) += 1;

      // Compute sizes of intersection withouth the observation obs
      counts_int(clust_minus[i], C_0_minus[i]) += 1;

    }
  }

  counts_set(find(counts_set > 0))= counts_set(find(counts_set > 0))/n;
  counts_int(find(counts_int > 0)) = counts_int(find(counts_int > 0))/n ;

  //--------------------------------------//

  for(int h = 0; h < H; h++){

    // Compute entropy of clustering
    H_c = -((counts_set(h) + 1/n)*log2(counts_set(h) + 1/n)) -
      sum(nonzeros(counts_set)%log2(nonzeros(counts_set))) +
      ((counts_set(h) > 0) ? (counts_set(h)*log2(counts_set(h))) : (0));


    // Compute mutual information
    H_cc = -(sum(nonzeros(counts_int)%log2(nonzeros(counts_int)))) +
      ((counts_int(h, C_0_minus(obs)) > 0) ? (counts_int(h, C_0_minus(obs))*log2(counts_int(h, C_0_minus(obs)))) : (0)) -
      ((counts_int(h, C_0_minus(obs)) +1/n)*log2(counts_int(h, C_0_minus(obs)) + 1/n));

    // VI distance (proportional up to H_0)
    dist(h) = 2*H_cc - H_c ;
  }

  return(dist);

}
