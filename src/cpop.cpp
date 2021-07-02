

#include "prune.impl.h"
#include "coeffs.update.h"

#include <Rcpp.h>
using namespace Rcpp;

#include <vector>

// [[Rcpp::plugins(cpp17)]]                                        

// [[Rcpp::export]]
std::vector<int> prune(const std::vector<double>& x,const int& nrows)
{
  std::vector<int> result(nrows);
  prune_impl(&x[0],&nrows,&result[0]);   
  return result;
}


// [[Rcpp::export]]
std::list<std::vector<double> > coeffs_update_cpp(const std::vector<double>& SXY,
						  const std::vector<double>& SX2,
						  const std::vector<double>& S,
						  const std::vector<double>& SS,
						  const std::vector<double>& Xs,
						  const std::vector<double>& SX,
						  const std::vector<double>& SP,
						  const std::vector<double>& Seglen,
						  const int& taustar,
						  const std::vector<int>& Sstar,
						  const double& Xt)  
{
  std::list<std::vector<double> > coeffs{calc_A(SX2,Xs,SX,SP,Seglen,taustar,Sstar),
					 calc_B(SX2,Xs,SX,SP,Seglen,taustar,Sstar,Xt),
					 calc_C(Xs,Seglen,taustar,Sstar,SXY,S),
					 calc_D(taustar,Sstar,SS),
					 calc_E(Xt,Seglen,taustar,Sstar,SXY,S),
					 calc_FF(SX2,SX,SP,Seglen,taustar,Sstar,Xt)};
  return coeffs;
}


