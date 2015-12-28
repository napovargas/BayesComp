#include <RcppArmadillo.h>
#include <cmath>
#include <time.h>
#include <updateThetaCpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

mat rwishart(int const& nu, mat const& V){

  int m = V.n_rows;
  mat T = zeros(m,m);
  mat R = zeros(m,m);
  for(int i = 0; i < m; i++) {
    T(i,i) = sqrt(rchisq(1,nu-i)[0]); 
  }
  
  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {    
      T(i,j) = rnorm(1)[0]; 
    }}
  
  mat C = trans(T)*chol(V);
  
    return R = trans(C)*C;
}

mat rinvwishart(int const& nu, mat const& V){

  int m = V.n_rows;
  mat T = zeros(m,m);
  mat IR = zeros(m,m);
  for(int i = 0; i < m; i++) {
    T(i,i) = sqrt(rchisq(1,nu-i)[0]);
  }
  
  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {    
      T(i,j) = rnorm(1)[0]; 
    }}
  
  mat C = trans(T)*chol(V);
  mat CI = solve(trimatu(C),eye(m,m)); 

    return IR = CI*trans(CI);
}

// [[Rcpp::export]] 
List BayesComp(mat M, mat Y, double nu, mat Omega, double gamma, uword nIter){
  Rcpp::Rcout << " ****************************************************** " << std::endl;
  Rcpp::Rcout << " Bayesian hierarchical model for unmixing compositions  " << std::endl;
  Rcpp::Rcout << " ****************************************************** " << std::endl;
  Rcpp::Rcout << "                                                        " << std::endl;
  Rcpp::Rcout << " Starting Gibbs sampler! " << std::endl;
  uword r             = M.n_rows;
  uword p             = M.n_cols;
  uword nAnim         = Y.n_cols;
  double ttime      = 0;
  clock_t start;
  clock_t end;
  vec init(p);
  mat Theta         = zeros(p, nAnim);//ones(p, nAnim)/p;
  for (uword t = 0; t < nAnim; t++){
   init          = runif(p);
   Theta.col(t) = init/accu(init);
  }
  mat store_Sigma   = zeros(nIter, r);
  mat store_Psi     = zeros(nIter, r);
  cube store_Theta  = zeros(nIter, p, nAnim);
  List Out;
  
  /* For updating Theta_i */
  vec ThetaStar;
  
  /* For updating Sigma */
  mat Sigma = rinvwishart(nAnim + nu, Omega);//eye(r, r);
  mat Psi   = rwishart(gamma, Omega);//eye(r, r);
  mat R     = zeros <mat> (r, r);
  
  start = clock();
  
  for(uword iter = 0; iter < nIter; iter++){
    for(uword i = 0; i < nAnim; i++){
      ThetaStar                                   = updateTheta(Theta.col(i), M, Sigma, Y.col(i));
      store_Theta(span(iter), span::all, span(i)) = ThetaStar;
      Theta.col(i)                                = ThetaStar;
      R                                           = R  + (Y.col(i) - M*ThetaStar)*trans(Y.col(i) - M*ThetaStar);
      }

      Sigma                 = rinvwishart(nAnim + nu, inv(R + Psi));
      Psi                   = rwishart(nu + gamma, Omega + Sigma);
      store_Sigma.row(iter) = trans(Sigma.diag()); 
      store_Psi.row(iter)   = trans(Psi.diag()); 
      R                     = zeros(r, r);
      if(iter % 1000 == 0){
        Rcpp::Rcout << " Iteration " << iter << " out of " << nIter << std::endl;
      }  
    }
    
    end = clock();
    ttime = ((double) (end - start)) / CLOCKS_PER_SEC;
    Rcpp::Rcout << nIter << " iterations in " << ttime << " seconds" << std::endl;
  
    Out["Theta"]        = store_Theta;
    Out["Sigma"]        = store_Sigma;
    Out["Psi"]          = store_Psi;
    Out["Elapsed time"] = ttime;
    return(wrap(Out));    
}