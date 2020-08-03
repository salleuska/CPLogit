#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Function which normalize probabilities in log scale
arma::vec normalize_probs(arma::vec x)
{
  arma::vec retval(x.n_elem);
  for(arma::uword u = 0; u < x.n_elem; ++u)
  {
    retval(u) =  exp(-log(sum(exp(x - x(u)))));
  }
  return retval;
}

// Function which relabel vector with contiguous entries from 0
// arma::uvec relabel(arma::uvec current){
//   arma::uvec labels = unique(current);
//   int n_groups = labels.n_elem;
//   arma::uvec index(n_groups);
//   arma::uvec res(current.n_elem, arma::fill::zeros);
//   arma::uvec tmp;
//   for(int j = 0; j < n_groups; ++j)
//   {
//     // For each of the labels find the first occurrence
//      tmp = find(current == labels(j),1);
//      index(j) = tmp(0);
//   }
//   // reored
//   index =  sort(index);
//
//   arma::uvec ind_relab;
//   // relabeling
//   for(int i = 0; i < index.n_elem; ++i){
//     ind_relab = find(current == current(index(i)));
//     res(ind_relab).fill(i);
//   }
//   return res;
// }

arma::uvec relabel(arma::uvec current)
{
  arma::uvec unique_vals = unique(current);
  int n_elem = unique_vals.n_elem;

  arma::uvec res(current.n_elem);
  arma::uvec index;
  for(int j = 0; j < n_elem; ++j)
  {
    index = find(current == unique_vals(j));
    res(index).fill(j);
  }

  return res;
}
