#include <RcppArmadillo.h>
#include <cmath>
#include <rpnorm1.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*// [[Rcpp::export]]*/
double rtnorm_MH(double X, double Mu, double Sigma, double Mum, double Mup){
  
  double Mu_new;
  double Mup_new;
  double Z;
  double Y;
  bool b0;
  Mu_new = Mu - Mum;
  Mup_new = Mup - Mum;
  if (Mu < Mup) {
    Z = rpnorm(Mu_new, Sigma);
    Z = Z + Mu_new;
  } else {
    Z = rpnorm(Mu_new, Sigma);
    Z = Z + -(Mu_new - Mup_new);
    Z = -(Z - Mup_new);
  }

  Z += Mum;
  if ((Z <= Mup) && (Z >= Mum)) {
    b0 = true;
  } else {
    b0 = false;
  }
  
  Y = Z * (double)b0 + X*(double)!b0;
  return(Y);
}
