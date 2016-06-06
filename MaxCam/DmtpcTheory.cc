#include "TMath.h"
#include "DmtpcTheory.hh"


double DmtpcTheory::ConvertXsectToRate( const Double_t xsect , const Double_t v_0, const Double_t A, const Double_t M_dark, const Double_t rho_D ) {

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


double DmtpcTheory::dRdERbetweenVearthandInfinity1DOverR0(Double_t *x, Double_t *par) {   
  // variables
   Double_t E_R=x[0];
  
  // user input parameters
  Double_t M_dark=par[0];
  Double_t M_target = par[2];
  Double_t v_0=par[3];
  Double_t v_earth=par[4];

  Double_t temp = 1/v_0;
  temp /= v_0;
  temp /= M_dark;
  temp /= (1./2.);
  // have to use temp, otherwise E_0 is too small
  // temp = 1/E_0

  // calculated quantities
  Double_t r = 4*M_dark*M_target/(TMath::Power((M_dark+M_target),2));
  Double_t v_min = TMath::Sqrt( E_R*temp/(r) )*v_0;
  
  Double_t value = temp;
  value *= 1/r;
  value *= sqrt(TMath::Pi())/4;
  value *= v_0/v_earth;
  value *= (TMath::Erf((v_min+v_earth)/v_0)-TMath::Erf((v_min-v_earth)/v_0));
  
  if( value< 0 ){
    value=0;
  }
  
  return value;
}

double  DmtpcTheory::dRdERbetweenVearthandInfinity1D(Double_t *x, Double_t *par) {
   

  // user input parameters
  Double_t M_dark=par[0];
  Double_t A=par[1];
  Double_t v_0=par[3];
  Double_t sigma_0=par[5];
  Double_t rho_D=par[6];
  x[0]=x[0]*1E-6; //convert keV to GeV
   
  Double_t value = ConvertXsectToRate( ( sigma_0*dRdERbetweenVearthandInfinity1DOverR0( x , par ) ) , v_0, A, M_dark, rho_D );
  
  x[0]=x[0]*1E6; //convert back to make ROOT play nice
  
  return value;
}

double DmtpcTheory::dRdERbetweenVearthandVesc1DOverR0(Double_t *x, Double_t *par) {

  // constants
   Double_t pi=TMath::Pi(); //pi
  
  // variables
   Double_t E_R=x[0];
  
  // user input parameters
  Double_t M_dark=par[0];
    Double_t M_target = par[2];
  Double_t v_0=par[3];
  Double_t v_earth = par[4];
  Double_t v_esc = par[7];

  Double_t temp = 1/v_0;
  temp /= v_0;
  temp /= M_dark;
  temp /= (1./2.);
  // have to use temp, otherwise E_0 is too small
  // temp = 1/E_0

  // calculated quantities
  Double_t r = 4*M_dark*M_target/(TMath::Power((M_dark+M_target),2));
  Double_t v_min = TMath::Sqrt( E_R*temp/(r) )*v_0;
  
  Double_t value = temp;
  value *= 1/r;
  value *= sqrt(TMath::Pi())/4;
  value *= v_0/v_earth;
  value *= (TMath::Erf((v_min+v_earth)/v_0)-TMath::Erf((v_min-v_earth)/v_0));
  if( value< 0 ){
    value=0;
  }

  double exponent = -v_esc*v_esc/(v_0*v_0);
  Double_t value2 = temp;
  value2 *=1/r;
  value2 *= exp(exponent);

  value -=value2;

  Double_t k1overk0 = TMath::Erf(v_esc/v_0)-2*v_esc/(sqrt(pi)*v_0)*exp(exponent);
  value /= k1overk0;

  if(value<0)
     value=0;
  
  return value;
}

double  DmtpcTheory::dRdERbetweenVearthandVesc1D(Double_t *x, Double_t *par) {
   

  // user input parameters
  Double_t M_dark=par[0];
  Double_t A=par[1];
  Double_t v_0=par[3];
  Double_t sigma_0=par[5];
  Double_t rho_D=par[6];
  x[0]=x[0]*1E-6; //convert keV to GeV
   
  Double_t value = ConvertXsectToRate( ( sigma_0*dRdERbetweenVearthandVesc1DOverR0( x , par ) ) , v_0, A, M_dark, rho_D );
  
  x[0]=x[0]*1E6; //convert back to make ROOT play nice
  
  return value;
}


double DmtpcTheory::dRdThOverR0(Double_t *x, Double_t *par) {
   
   // variables
   Double_t E_R=par[7]*1E-6;
   Double_t cosPsi=x[0];
   
   // user input parameters
   Double_t M_dark=par[0];
   Double_t M_target = par[2];
   Double_t v_0=par[3];
   Double_t v_earth =par[4];
   
   Double_t temp = 1/v_0;
   temp /= v_0;
   temp /= M_dark;
   temp /= (1./2.);
   // have to use temp, otherwise E_0 is too small
   
   // calculated quantities
   Double_t r = 4*M_dark*M_target/(TMath::Power((M_dark+M_target),2));
   Double_t v_min = TMath::Sqrt( E_R*temp/(r) )*v_0;
   
   Double_t exponent=( v_earth*cosPsi - v_min );
   exponent *= ( v_earth*cosPsi - v_min );
   exponent /= v_0;
   exponent /= v_0;
   
   Double_t value = TMath::Exp( -exponent );
   value *= ( 1. )/( 2.);
   value *= temp;
   value *= 1/r;
   
   if( value< 0 ){
      value=0;
   }
   
   return value;
}

double DmtpcTheory::dRdTh(Double_t *x, Double_t *par) {
   
  // user input parameters
  Double_t M_dark=par[0];
  Double_t A=par[1];
  Double_t v_0=par[3];
  Double_t sigma_0=par[5];
  Double_t rho_D=par[6];
    
  Double_t value = ConvertXsectToRate( ( sigma_0*dRdThOverR0( x , par ) ) , v_0, A, M_dark, rho_D );
  
   return value;
}

double DmtpcTheory::dRdERbetweenVearthandInfinity2DOverR0(Double_t *x, Double_t *par) {
  
  // variables
   Double_t E_R=x[0];
   Double_t cosPsi=x[1];
   
  // user input parameters
  Double_t M_dark=par[0];
  Double_t M_target=par[2];
  Double_t v_0=par[3];
  Double_t v_earth = par[4];

  Double_t temp = 1/v_0;
  temp /= v_0;
  temp /= M_dark;
  temp /= (1./2.);
  // have to use temp, otherwise E_0 is too small

  // calculated quantities
  Double_t r = 4*M_dark*M_target/(TMath::Power((M_dark+M_target),2));
  Double_t v_min = TMath::Sqrt( E_R*temp/(r) )*v_0;
  
  Double_t exponent=( v_earth*cosPsi - v_min );
  exponent *= ( v_earth*cosPsi - v_min );
  exponent /= v_0;
  exponent /= v_0;
  
  Double_t value = TMath::Exp( -exponent );
  value *= ( 1. )/( 2.);
  value *= temp;
  value *= 1/r;
  
  if( value< 0 ){
    value=0;
  }
  
  return value;
}




double  DmtpcTheory::dRdERbetweenVearthandInfinity2D(Double_t *x, Double_t *par) {
   
  // user input parameters
  Double_t M_dark=par[0];
  Double_t A=par[1];
  Double_t v_0=par[3];
  Double_t sigma_0=par[5];
  Double_t rho_D=par[6];
  x[0]=x[0]*1E-6; //convert keV to GeV
   
  Double_t value = ConvertXsectToRate( ( sigma_0*dRdERbetweenVearthandInfinity2DOverR0( x , par ) ) , v_0, A, M_dark, rho_D );
  
  x[0]=x[0]*1E6; //convert back to make ROOT play nice
  
  return value;
}
