#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

/*// [[Rcpp::export]]*/

arma::uvec sequence(arma::uword from, arma::uword to){
    arma::uword nelems = to - from + 1;
    arma::uvec x(nelems);
    for (arma::uword i = 0; i < nelems; i++){
        x(i) = from + i;
    }
    return x;
}