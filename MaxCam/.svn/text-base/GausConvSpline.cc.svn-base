#include "GausConvSpline.hh" 
#include "TMath.h"
#include <cmath>
#include <stdio.h>

static const int debug= 0; 


GausConvSpline::GausConvSpline(double s, TSpline3 * spline, double l, double u) 
{

  lower.push_back(l); 

  for (int i = 0; i < spline->GetNp()-1; i++)
  {
    double xlow, xhigh; 
    double yignore; 
    spline->GetKnot(i,xlow,yignore); 
    spline->GetKnot(i+1,xhigh,yignore); 

    double bound = xlow + xhigh; 
    bound/=2; 
    if (bound > l && bound < u)
    {
      upper.push_back(bound); 
      lower.push_back(bound); 
    }
    else if (bound >= u) 
      break; 
  }

  upper.push_back(u); 

  int knot = spline->FindX(l);  

  double Xa,Ya,a,b,c,A,B; 
  double sqrt2pi = sqrt(2*M_PI); 
  for (unsigned i = 0; i < lower.size(); i++)
  {
    
    coeffs_lower_gaus.push_back(vector<double>(3,0)); 
    coeffs_upper_gaus.push_back(vector<double>(3,0)); 
    coeffs_erf.push_back(vector<double>(4,0)); 

    spline->GetCoeff(knot,Xa,Ya,a,b,c); 
    A = lower[i]; 
    B = upper[i]; 
    
    if (debug) printf("Working on interval (%f,%f)\n",A,B); 
    if (debug)  printf("\t (Xa,Ya,a,b,c) = (%f,%f,%f,%f,%f)\n", Xa,Ya,a,b,c); 

    coeffs_lower_gaus[i][0] = c*s/sqrt2pi; 
    coeffs_lower_gaus[i][1] = s/sqrt2pi * ( b + A*c - 3*c*Xa); 
    coeffs_lower_gaus[i][2] = s/sqrt2pi * (a + A*b + A*A*c + 2*c*s*s - 2*b*Xa - 3*A*c*Xa + 3*c*Xa*Xa); 

    coeffs_upper_gaus[i][0] = c*s/sqrt2pi; 
    coeffs_upper_gaus[i][1] = s/sqrt2pi * ( b + B*c - 3*c*Xa); 
    coeffs_upper_gaus[i][2] = s/sqrt2pi * (a + B*b + B*B*c + 2*c*s*s - 2*b*Xa - 3*B*c*Xa + 3*c*Xa*Xa); 

    coeffs_erf[i][0] = c/2; 
    coeffs_erf[i][1] = 1./2 * (b - 3*c*Xa); 
    coeffs_erf[i][2] = 1./2 * (a + 3*c*s*s - 2*b*Xa + 3*c*Xa*Xa); 
    coeffs_erf[i][3] = 1./2 * (b*s*s - a*Xa - 3*c*s*s*Xa + b*Xa*Xa - c * Xa*Xa*Xa + Ya); 

    knot++;
  }

  sgma = s; 
  sgma2 = 2*s*s; 
  sqrt2s = sqrt(2) * s; 

}

double GausConvSpline::AntiDerivative( double *x, double *p)
{
  double T = *x; 

}

double GausConvSpline::Eval (double * x, double * p) 
{
  double T = *x; 
  double lb = T - *p * sgma; 
  double ub = T + *p * sgma; 

  if (debug) printf("T = %f\n, bounds = (%f,%f)",T,lb,ub); 

  if (ub < lower[0]) return 0; 
  if (lb > upper[upper.size()-1]) return 0; 

  int first,last; 
  double T2 = T*T; 
  double T3 = T*T*T; 

  if (lb < lower[0]) first = 0; 
  else
  {
    first = TMath::BinarySearch(lower.begin(), lower.end(), lb) - lower.begin(); 
  }

  last = TMath::BinarySearch(upper.begin(), upper.end(), ub) - upper.begin() + 1; 
  if (last == (int) upper.size()) last--; 

  if (debug) printf("Using intervals %d to %d\n",first,last); 

  double answer = 0; 

  for (int i = first; i <= last; i++)
  {
    double lower_gaus = exp(-pow(lower[i] - T,2)/sgma2);   
    double upper_gaus = exp(-pow(upper[i] - T,2)/sgma2);   
    double erf = (TMath::Erf((T-lower[i])/sqrt2s) - TMath::Erf((T-upper[i])/sqrt2s)); 
    if (debug) printf("Interval %d, lower_gaus, upper_gaus, erf: %f,%f,%f\n", i, lower_gaus, upper_gaus, erf); 
    answer += lower_gaus * (T2 * coeffs_lower_gaus[i][0] + T * coeffs_lower_gaus[i][1] + coeffs_lower_gaus[i][2]); 
    answer -= upper_gaus * (T2 * coeffs_upper_gaus[i][0] + T * coeffs_upper_gaus[i][1] + coeffs_upper_gaus[i][2]); 
    answer += erf * ( T3 * coeffs_erf[i][0] + T2 * coeffs_erf[i][1] + T * coeffs_erf[i][2] + coeffs_erf[i][3]); 
  }

  if (debug) printf("Answer: %f\n\n", answer); 
  return answer; 

}
