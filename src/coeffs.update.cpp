
#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <list>
#include <algorithm>

// [[Rcpp::plugins(cpp17)]]                                        
std::vector<double> calc_A(const std::vector<double>& SX2,
			   const std::vector<double>& Xs,
			   const std::vector<double>& SX,
			   const std::vector<double>& SP,
			   const std::vector<double>& Seglen,
			   const int& taustar,
			   const std::vector<int>& Sstar)
{
  std::vector<double> A(Sstar.size());
  std::transform(Xs.begin(),Xs.end(),Sstar.begin(),A.begin(),

		 [&SX2,&SX,&Xs,&taustar,&SP](auto& xs, auto& sstar)
		 {
		   return SX2[taustar] - SX2[sstar] - 2*xs*(SX[taustar] - SX[sstar]) + (SP[taustar] - SP[sstar])*xs*xs;
		 });
  std::transform(Seglen.begin(),Seglen.end(),A.begin(),A.begin(),
		 [](auto& s,auto& a)
		 {
		   return a/(s*s);
		 });
  return A;
}


std::vector<double> calc_B(const std::vector<double>& SX2,
			   const std::vector<double>& Xs,
			   const std::vector<double>& SX,
			   const std::vector<double>& SP,  
			   const std::vector<double>& Seglen,
			   const int& taustar,
			   const std::vector<int>& Sstar,
			   const double& Xt)
{
  std::vector<double> B(Sstar.size());

  std::transform(Xs.begin(),Xs.end(),Sstar.begin(),B.begin(),

		 [&SX2,&SX,&Xs,&taustar,&SP,&Xt](auto& xs, auto& sstar)
		 {

		   return (Xt+xs)*(SX[taustar] - SX[sstar]) - (SP[taustar] - SP[sstar])*Xt*xs - (SX2[taustar] - SX2[sstar]);
		 });
  std::transform(Seglen.begin(),Seglen.end(),B.begin(),B.begin(),
		 [](auto& s,auto& b)
		 {
		   return 2*b/(s*s);
		 });      
  return B;
  
}


std::vector<double> calc_C(const std::vector<double>& Xs,
			   const std::vector<double>& Seglen,
			   const int& taustar,
			   const std::vector<int>& Sstar,
			   const std::vector<double>& SXY,
			   const std::vector<double>& S)
{
  std::vector<double> C(Sstar.size());
  std::transform(Xs.begin(),Xs.end(),Sstar.begin(),C.begin(),

		 [&SXY,&S,&Xs,&taustar](auto& xs, auto& sstar)
		 {
		   return SXY[taustar] - SXY[sstar] - xs*(S[taustar] - S[sstar]);
		 });
    std::transform(Seglen.begin(),Seglen.end(),C.begin(),C.begin(),
		 [](auto& s,auto& c)
		 {
		   return -(2.0/s)*c;
		 });      

  return C;
}


std::vector<double> calc_D(const int& taustar,
			   const std::vector<int>& Sstar,
			   const std::vector<double>& SS)
{
  std::vector<double> D(Sstar.size());
  std::transform(Sstar.begin(),Sstar.end(),D.begin(),
		 [&SS,&taustar](auto& sstar)
		 {
		   return SS[taustar] - SS[sstar];
		 });
  return D;
  
}

std::vector<double> calc_E(const double& Xt,
			   const std::vector<double>& Seglen,
			   const int& taustar,
			   const std::vector<int>& Sstar,
			   const std::vector<double>& SXY,
			   const std::vector<double>& S)
{
  std::vector<double> E(Sstar.size());
  std::transform(Sstar.begin(),Sstar.end(),E.begin(),
		 [&SXY,&S,&Xt,&taustar](auto& sstar)
		 {
		   return Xt*(S[taustar] - S[sstar]) - (SXY[taustar] - SXY[sstar]);
		 });
    std::transform(Seglen.begin(),Seglen.end(),E.begin(),E.begin(),
		 [](auto& s,auto& e)
		 {
		   return -(2.0/s)*e;
		 });      

  return E;
}

std::vector<double> calc_FF(const std::vector<double>& SX2,
			   const std::vector<double>& SX,
			   const std::vector<double>& SP,
			   const std::vector<double>& Seglen,
			   const int& taustar,
			    const std::vector<int>& Sstar,
			    const double& Xt)
{
  std::vector<double> FF(Sstar.size());
  std::transform(Sstar.begin(),Sstar.end(),FF.begin(),

		 [&SX2,&SX,&Xt,&taustar,&SP](auto& sstar)
		 {
		   return SX2[taustar] - SX2[sstar] - 2*Xt*(SX[taustar] - SX[sstar]) + (SP[taustar] - SP[sstar])*Xt*Xt;
		 });
  std::transform(Seglen.begin(),Seglen.end(),FF.begin(),FF.begin(),
		 [](auto& s,auto& ff)
		 {
		   return ff/(s*s);
		 });
  return FF;
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










