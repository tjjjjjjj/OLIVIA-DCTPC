#ifndef GAUS_CONV_SPLINE_HH
#define GAUS_CONV_SPLINE_HH

#include "TSpline.h" 
#include "TObject.h" 
#include <vector> 

using namespace std; 

class GausConvSpline : public TObject
{
  public:
    GausConvSpline(double sigma, TSpline3 * spline, double a, double b); 

    // 1 dim, takes one parameter (how many sigma away to consider segments)
    double Eval(double * x, double * p); 
    double AntiDerivative(double * x, double * p); 

  private: 
    vector<double> lower;  //lower bounds
    vector<double> upper;  //uppper bounds 
    vector<vector<double> > coeffs_lower_gaus; 
    vector<vector<double> > coeffs_upper_gaus; 
    vector<vector<double> > coeffs_erf; 
    double sgma; 
    double sgma2; 
    double sqrt2s; 


};


#endif


