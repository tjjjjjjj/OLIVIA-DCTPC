#include "SimChamber.hh"
#include "TVector.h"
#include "TVectorD.h"
#include "TGraph.h"
#include <iostream>
//__________________________________________
/* Begin_Html
<center><h2>The SimChamber Class </h2></center>
This class is essentially just a basic C++ structure with set and get commands.
The only other functions of the class are the constructor and a function 
to get the z-dependence of the diffusion using data taken from the diffusion
study paper. 
End_Html*/
//__________________________________________

ClassImp(SimChamber)
SimChamber::SimChamber()
{
  //>>>>>Default Constructor <<<<<
  fHeight = 20*25.4;//20 inch
  fDriftLength = 197;//19.7 cm
  fDriftVoltage = 5000;
  fAnodeVoltage = 720; 
  fTemperature = 300;
  fPressure = 75;
  fDOverMu = 0.00419*fDriftVoltage / (2*fDriftLength);
  fDiffusionConstTerm = 0.25;
  fElecScintillation = 1.0;
  fNuclScintillation = 0.3;
  fAttenuation = 0;
  fElectronPerkeV = 10000;
  fLongDiffConstTerm = 400;
  fLongDiffDzTerm = 3;
  fOrientationAngle = 0;
}

void
SimChamber::setDOverMu()
{
  //Set the value of D / mu using the data from the diffusion study paper
  //Depends only on drift voltage, drift length, temperature and pressure
  double kb = 1.38065e-19;//Boltzmann in kg-cm^2/s^2
  double toTd = 1e17;//Convert from V-cm^2 to Td
  double torr = 0.133322;//1 torr in kg mm^-1 s^-2
  TVectorD E_N(14);
  TVectorD D_mu(14);
  for (int i = 0; i < 14; i++) E_N[i] = i+1;
  D_mu[0] = 0.027;
  D_mu[1] = 0.027;
  D_mu[2] = 0.032;
  D_mu[3] = 0.033;
  D_mu[4] = 0.033;
  D_mu[5] = 0.039;
  D_mu[6] = 0.037;
  D_mu[7] = 0.042;
  D_mu[8] = 0.044;
  D_mu[9] = 0.051;
  D_mu[10] = 0.055;
  D_mu[11] = 0.055;
  D_mu[12] = 0.062;
  D_mu[13] = 0.078;
  TGraph gr(E_N,D_mu);
  //E_N = Electric field / number density = V/L * kT / P for ideal gas
  //BE CAREFUL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //double e_n = fDriftVoltage * fTemperature / (fDriftLength * fPressure) * kb/torr * toTd;

  double e_n = fDriftVoltage * fTemperature / (fDriftLength * 1140.) * kb/torr * toTd;
  std::cout<<fDriftVoltage<<" "<<fTemperature<<" "<<fDriftLength<<" "<<fPressure<<" "<<e_n<<" "<<gr.Eval(e_n)<<std::endl;

  fDOverMu = gr.Eval(e_n);
}
