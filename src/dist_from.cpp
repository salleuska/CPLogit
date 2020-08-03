#include <RcppArmadillo.h>
// Using code from fxt library
// CREDITS TO GIVE
#include  "dist_from.h"

// Function which compute the VI distance from c0
// and all the other partition in the namespace
// Generate all the space by using code form fxt library

// Return also blocks sizes for each partition (suff stats for eppf)
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// double dp_eppf(arma::vec clust, double alpha){
//   // Output
//   double log_prob;
//
//   int n = clust.n_elem;
//   arma::uvec clust_labels= unique(clust);
//   int n_blocks = clust_labels.n_elem;
//
//   double sum_lgamma_blocks = 0;
//   // Block sizes
//   for(int h = 0; h < n_blocks; ++h){
//     arma::uvec sel = find(clust == h);
//     sum_lgamma_blocks += lgamma(sel.n_elem);
//   }
//
//   log_prob = n_blocks*log(alpha) + lgamma(alpha)  -lgamma(alpha + n) + sum_lgamma_blocks;
//
//   return exp(log_prob);
// }

// Bell numbers - 0 to 25 - sequence A000110
unsigned long long bell_array[] = { 1,1,2,5,15,52,203,877,4140,21147,115975,678570, 4213597,27644437,190899322,1382958545 ,10480142147, 82864869804,682076806159,5832742205057, 51724158235372,474869816156751,4506715738447323, 44152005855084346,445958869294805289, 4638590332229999353};

// [[Rcpp::export]]
List dist_from(arma::uvec clust, bool return_partitions = false)
{
  // NEED TO ADD CHECKS FOR N WHEN DOING A LIBRARY
  // number of observations
  int n = clust.n_elem;
  int n_tot =  bell_array[n];
  arma::vec out(n_tot, arma::fill::zeros);
  arma::mat blocks(n_tot, n, arma::fill::zeros);
  arma::umat partitions(n_tot, n, arma::fill::zeros);
  arma::uvec c(n, arma::fill::zeros);
  // object with all set partitions of n
  setpart P(n);

  for(int i = 0; i < n_tot; ++i)
  {
    for (int j = 1; j < n +1; ++j){
      c(j-1) = P.as_[j];
    }
    out(i) = vi_distC(c, clust);
    // Fills matrix with partitions
    partitions.row(i) = c.t();

    // Fills matrix with block sizes
    arma::uvec clust_labels= unique(clust);
    int n_blocks = clust_labels.n_elem;
    for(int h = 0; h < n; ++h){
      arma::uvec sel = find(c == h);
      blocks(i, h) = sel.n_elem;
    }
    P.next();
  }

  if(return_partitions == true) {
    return List::create(_["distances"] = out,
                      _["block_sizes"] = blocks,
                      _["partitions"] = partitions);
  }
  return List::create(_["distances"] = out,
                    _["block_sizes"] = blocks);
}
