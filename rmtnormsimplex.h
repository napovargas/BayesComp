#include <RcppArmadillo.h>
#include <cmath> 
#include <sequence.h>
#include <rtnormcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

vec rmtnorm_simplex(vec S, vec Mean, mat Var){
  vec Mu            = Mean;
  uword p           = Mu.n_elem;
  uword j           = 0;
  vec W             = zeros(p);
  vec Mu_sv         = zeros(p);
  vec Var_sv        = zeros(p);
  vec Sd_sv         = zeros(p);
  vec Sj;
  vec Muj;
  vec Rv(p - 1);
  uvec P            = sequence(1, p);
  uvec shuffledP    = shuffle(P);
  mat vecSigma      = zeros(p - 1, p);
  mat Rm            = zeros(p, p);
  cube matSigma     = zeros(p - 1, p - 1, p);
  
  for(uword r = 0;r < p; r++){
    Rm = Var;
    Rm.shed_row(r);
    Rv = Rm.col(r);
    Rm.shed_col(r);
    matSigma.slice(r) = inv(Rm);
    vecSigma.col(r) = Rv;
  } 
    for (uword iter = 0; iter < 10; iter ++){
        for(uvec::iterator jit = shuffledP.begin(); jit != shuffledP.end(); jit++){
            j = (int)*jit - 1;
            Sj = S;
            Sj.shed_row(j);
            Muj = Mu;
            Muj.shed_row(j);
            Mu_sv(j) = (Mu(j) + (trans(vecSigma.col(j))*matSigma.slice(j)*(Sj - Muj)))[0];
            Var_sv(j) = (Var(j, j) - (trans(vecSigma.col(j))*matSigma.slice(j)*vecSigma.col(j)))[0];
            Sd_sv(j) = sqrt(std::abs(Var_sv(j)));
            W(j) = rtnorm(0.0, std::min(1.0, std::abs(1 - sum(S) + S(j))), Mu_sv(j), Sd_sv(j));
            //W(j) = rtnorm_MH(S(j), Mu_sv(j), Sd_sv(j), 0.0, (1 - sum(S) + S(j)));

        }
    }
  return(W);
}
