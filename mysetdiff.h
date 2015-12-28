#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

/*// [[Rcpp::export]]*/
arma::uvec mysetdiff(arma::uvec& x, arma::uvec& y){

  x = arma::unique(x);
  y = arma::unique(y);
  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uvec q1 = arma::find(x == y[j]);
    x.shed_row(q1(0));
  }

  return x;
}