#include "TMath.h"

Double_t ConvertXsectToRate( const Double_t xsect , const Double_t v_0, const Double_t A, const Double_t M_dark, const Double_t rho_D ) {

  const Double_t N0_over_1amu=3.62661882e50; // kg^-1
  const Double_t day_over_sec=86400;
  const Double_t km_over_cm=100000;
  const Double_t pi=TMath::Pi();
  const Double_t c = 3.e5; // in km/s

  Double_t rate = xsect;
  rate *= ( 2.)/( TMath::Sqrt(pi) );
  rate /= A;
  rate /= M_dark/rho_D;
  rate *= ( v_0*c );
  rate *= N0_over_1amu;
  rate /= day_over_sec;
  rate /= km_over_cm;
  
  return rate;
}

Double_t ConvertRateToXsect( const Double_t rate , const Double_t v_0, const Double_t A, const Double_t M_dark, const Double_t rho_D ) {
  
  const Double_t N0_over_1amu=3.62661882e50; // kg^-1
  const Double_t day_over_sec=86400;
  const Double_t km_over_cm=100000;
  const Double_t pi=TMath::Pi();
  const Double_t c = 3.e5; // in km/s

  Double_t xsect = rate;
  xsect *= TMath::Sqrt(pi)/2.;
  xsect *= A;
  xsect *= M_dark/rho_D;
  xsect /= ( v_0*c );
  xsect /= N0_over_1amu;
  xsect *= day_over_sec;
  xsect *= km_over_cm;
  
  return xsect;
} 

