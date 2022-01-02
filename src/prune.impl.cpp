

#include <math.h>
#include <stdio.h>
#include <limits.h>
#include "R.h"


#include "Rcpp.h"
using namespace Rcpp;


#include <cmath>
#include <limits>



bool cmp_double(const double& a, const double& b)
{
  return a == b;
  /*
  static constexpr auto epsilon = std::numeric_limits<double>::epsilon();
    if (std::abs(a - b) <= epsilon)
        return true;
    return std::abs(a - b) <= epsilon * std::max(std::abs(a), std::abs(b));
  */
}


// need to put const n this
void prune_impl(const double* x, const int* nrows, int* Sets)
// void prune2R_orig_impl(const *x, int *nrows, int *Sets) ---- original
{
  
  int i;
  double tcurr= -DBL_MAX;
  int whichfun=0;
  double minA=*(x+2), maxB=*(x+3), minC=*(x+4);
    for ( i = 1 ; i < *nrows ; i++ ) 
    {
      if ( *(x+5*i+2) < minA ) 
      {
        minA = *(x+5*i+2);
	maxB = *(x+5*i+3);
	minC = *(x+5*i+4);
        whichfun = i;
      }
      // else if(*(x+5*i+2) == minA){
      else if(cmp_double(*(x+5*i+2),minA)){
	if(*(x+5*i+3) > maxB){
	  maxB = *(x+5*i+3);
	  minC = *(x+5*i+4);
	  whichfun = i;
	}
	// else if(*(x+5*i+3) == maxB){
	else if(cmp_double(*(x+5*i+3),maxB)){
	  if(*(x+5*i+4) < minC){
	    minC = *(x+5*i+4);
	    whichfun = i;
	  }
	}
      }
    }

 std::vector<int> logicint(*nrows);  
 // int logicint[*nrows]; /*1 indicates NA*/
 for(i=0;i<*nrows;i++){
   logicint[i]=1;
   *(Sets+i)=0;
 }
 *(Sets+whichfun)=1;
 int sum=*nrows;
    while(sum>1){
      std::vector<double> intercepts(*nrows,0.0);  
      // double intercepts[*nrows];
      logicint[whichfun]=0;
      intercepts[whichfun]=0;
      for (i=0;i<*nrows;i++){
	if(logicint[i]!=0){
	  double A=*(x+5*i+2)-*(x+5*whichfun+2);
	  double B=*(x+5*i+3)-*(x+5*whichfun+3);
	  double C=*(x+5*i+4)-*(x+5*whichfun+4); /*creates the diffcoeff function*/
	  double disc=B*B-4*A*C;
	  if(disc<0){
	    intercepts[i]=0/*NA*/;
	    logicint[i]=0;
	  }
	  else{
	    // if(A==0){
	    if(cmp_double(A,0)){
	      // if(B==0){
	      if(cmp_double(B,0)){
		intercepts[i]=0;
		logicint[i]=0;
		}
	      else if((-C/B)>tcurr){
		intercepts[i]=-C/B;
		logicint[i]=1;
	      }
	      else{
		intercepts[i]=0/*NA*/;
		logicint[i]=0;
	      }
	    }
	    else{
	      double vec2a=(-B-sqrt(disc))/(2*A);
	      double vec2b=(-B+sqrt(disc))/(2*A);
	      if(vec2a<=tcurr && vec2b<=tcurr){
		intercepts[i]=0/*NA*/;
		logicint[i]=0;
	      }
	      else if(vec2a<=tcurr && vec2b>tcurr){
		intercepts[i]=vec2b;
		logicint[i]=1;
	      }
	      else if((vec2a>tcurr && vec2b<=tcurr) || (vec2a>tcurr && vec2b>tcurr && vec2a<vec2b)){
		intercepts[i]=vec2a;
		logicint[i]=1;
	      }
	      else{
		intercepts[i]=vec2b;
		logicint[i]=1;
	      }
	    }
	  }
	}
      }
      int j;
      sum=0;
      int whichfunnew = 0;
      // double minimum=LONG_MAX;
      double minimum=DBL_MAX;
      for(j=0;j<*nrows;j++){
	if(logicint[j]!=0){
	  sum=sum+logicint[j];
	  if ( intercepts[j] < minimum ) {
	    minimum = intercepts[j];
	    whichfunnew = j;
	  }
	}
      }
      logicint[whichfun]=1;
      if(sum!=0){
	tcurr = minimum;
	whichfun = whichfunnew;
	*(Sets+whichfun)=1;
      }
    }
}





