#ifndef ___COEFFS_UPDATE_H___
#define ___COEFFS_UPDATE_H___

#include <vector>


std::vector<double> calc_A(const std::vector<double>&,
			   const std::vector<double>&,
			   const std::vector<double>&,
			   const std::vector<double>&,
			   const std::vector<double>&,
			   const int&,
			   const std::vector<int>&);



std::vector<double> calc_B(const std::vector<double>&,
			   const std::vector<double>&,
			   const std::vector<double>&,
			   const std::vector<double>&,  
			   const std::vector<double>&,
			   const int&,
			   const std::vector<int>&,
			   const double&);



std::vector<double> calc_C(const std::vector<double>&,
			   const std::vector<double>&,
			   const int&,
			   const std::vector<int>&,
			   const std::vector<double>&,
			   const std::vector<double>&);



std::vector<double> calc_D(const int&,
			   const std::vector<int>&,
			   const std::vector<double>&);


std::vector<double> calc_E(const double&,
			   const std::vector<double>&,
			   const int&,
			   const std::vector<int>&,
			   const std::vector<double>&,
			   const std::vector<double>&);


std::vector<double> calc_FF(const std::vector<double>&,
			   const std::vector<double>&,
			   const std::vector<double>&,
			   const std::vector<double>&,
			   const int&,
			    const std::vector<int>&,
			    const double&);

#endif
