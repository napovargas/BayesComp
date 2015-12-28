#include <RcppArmadillo.h>
#include <cmath>
#include <mysetdiff.h>
#include <rmtnormsimplex.h>
//#include <rtnormcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]] 
arma::vec updateTheta(arma::vec ThetaStar, arma::mat M, arma::mat Sigma, arma::vec y){
    uword p                 = M.n_cols;
    uvec plants;
    uvec tmp;
    uvec lastj;
    uvec inTheta;
    uvec estim;
    uvec last;
    vec Theta_p;
    mat Mstar;
    mat mp;
    mat Var;
    vec Mean                = zeros(p - 1);
    vec Mu                  = zeros(p - 1);
    mat One                 = ones <mat> (p - 1, 1);
    plants                  = sequence(1, p);
    tmp                     = shuffle(plants);
    lastj                   = tmp(0);
    inTheta                 = mysetdiff(tmp, lastj);
    estim                   = inTheta -  1;
    last                    = lastj - 1;
    Mstar                   = M.cols(estim);
    mp                      = M.cols(last);
    Var                     = (trans(Mstar - mp*trans(One))*solve(Sigma, Mstar - mp*trans(One)));
    if (p == 2){
        Mean                    = ((1/Var[0,0])*trans(Mstar - mp*trans(One))*solve(Sigma, y - mp));
        Mu                      = rtnorm(0.0, 1.0, Mean(0), sqrt(1/Var[0, 0]));
        ThetaStar.elem(estim)   = Mu;
        ThetaStar.elem(last)    = 1 - Mu;
    }
    else {
        Mean                    = (inv(Var)*trans(Mstar - mp*trans(One))*solve(Sigma, y - mp));
        Mu                      = rmtnorm_simplex(ThetaStar.elem(estim), Mean, solve(Var, eye(p - 1, p - 1)));
        Theta_p << 1 - accu(Mu);
        ThetaStar.elem(estim)   = Mu;
        ThetaStar.elem(last)    = Theta_p;
    }
    /*Rcpp::Rcout << "Mean: " << std::endl;
    Rcpp::Rcout << Mean << std::endl;
    Rcpp::Rcout << "Mu: " << std::endl;
    Rcpp::Rcout << Mu << std::endl;
    Rcpp::Rcout << "Var: " << std::endl;
    Rcpp::Rcout << Var << std::endl;*/

    return (ThetaStar);
}