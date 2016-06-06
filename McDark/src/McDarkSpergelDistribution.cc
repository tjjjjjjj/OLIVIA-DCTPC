#include "TMath.h"

#include "McDarkSpergelDistribution.hh"
#include "differentialRateSpergel.hh"
#include "Randomize.hh"
#include "globals.hh"


McDarkSpergelDistribution::McDarkSpergelDistribution()
{
  // various variables to configure plotting
  // the differential rate will be computed between these two limits of nuclear kinetic recoil energy, in GeV.  
  const Double_t lowerERLimit=0./1.e6; // to put into GeV
  const Double_t upperERLimit=25./1.e6; // to put into GeV
  const Double_t lowerCosPsiLimit=-1.0;
  const Double_t upperCosPsiLimit=1.0;

  differentialRateSpergel = new TF2("dRdERbetweenVearthandInfinity2DERProjection",
                                    dRdERbetweenVearthandInfinity2D,lowerERLimit,upperERLimit,lowerCosPsiLimit,upperCosPsiLimit,10);

  // the velocity of light m/s
  const Double_t c = 3.e8;
  // WIMP mass in GeV
  Double_t M_dark=60.;
  // atomic number of the target.  132=Xenon
  Double_t A=132.;
  // the average dark matter velocity (unitless)
  // unless I say otherwise, these values are taken from the Lewin and Smith paper
  Double_t v_0=230.e3/c;
  // the density of WIMPs in the galaxy halo
  Double_t rho_D=0.4;  // GeV per cm^3 c^2
  // the crossection of dark matter interaction with the stuff in detector.  Can be anything
  Double_t sigma_0=1e-44; // cm^2

  // the escape velocity of dark matter in the galaxy.  If dark matter has a local velocity greater than this
  // it boils out of the galaxy.
  //Double_t v_esc=600.e3/c; // comment to avoid warning [DD 071011]

  // the velocity of the Earth around the sun
  Double_t v_E_0=244.e3/c;
  // the variance of the Earth's velocity around the sun, depending on where the Earth is around the sun
  // the total earth velocity looks like v_E=v_E_0+delta_v_E_0*Sin(y*2*PI)
  Double_t delta_v_E_0=15.e3/c;
  // per above comment, y indexes the time of year; i.e., where the Earth is around the sun.  See the
  // correspondences above.  I'll just set it to zero
  Double_t y=0;

  // set parameters
  differentialRateSpergel->SetParameters(M_dark,A,v_0,v_E_0,delta_v_E_0,y,sigma_0,rho_D);
}




void McDarkSpergelDistribution::GetRandomValues(Double_t e_and_psi[]) {

    TF2 *fun=differentialRateSpergel;
    
   //  Check if integral array must be build
   Int_t i,j,cell;
   Double_t dx   = (fun->GetXmax()-fun->GetXmin())/fun->GetNpx();
   Double_t dy   = (fun->GetYmax()-fun->GetYmin())/fun->GetNpy();
   Int_t ncells = fun->GetNpx()*fun->GetNpy();
   Double_t *fIntegral=0;
   if (fIntegral == 0) {
      fIntegral = new Double_t[ncells+1];
      fIntegral[0] = 0;
      Double_t integ;
      Int_t intNegative = 0;
      cell = 0;
      for (j=0;j<fun->GetNpy();j++) {
          for (i=0;i<fun->GetNpx();i++) {
              integ = fun->Integral(fun->GetXmin()+i*dx,fun->GetXmin()+i*dx+dx,fun->GetYmin()+j*dy,fun->GetYmin()+j*dy+dy);
              if (integ < 0) {intNegative++; integ = -integ;}
              fIntegral[cell+1] = fIntegral[cell] + integ;
              cell++;
          }
      }
      if (intNegative > 0) {
          G4cerr << "McDarkSpergelDistribution: probability function has negative values - abs assumed";
      }
      if (fIntegral[ncells] == 0) {
          G4cerr << "McDarkSpergelDistribution: Integral of function is zero" << G4endl;
          return;
      }
      for (i=1;i<=ncells;i++) {  // normalize integral to 1
          fIntegral[i] /= fIntegral[ncells];
      }
   }

   
// return random numbers
   Double_t r,ddx,ddy,dxint;
   r     = G4UniformRand();
   cell  = TMath::BinarySearch(ncells,fIntegral,r);
   dxint = fIntegral[cell+1] - fIntegral[cell];
   if (dxint > 0) ddx = dx*(r - fIntegral[cell])/dxint;
   else           ddx = 0;
   ddy = dy*G4UniformRand();
   j   = cell/fun->GetNpx();
   i   = cell%fun->GetNpx();
   e_and_psi[0] = fun->GetXmin() +dx*i +ddx;
   e_and_psi[1] = fun->GetYmin() +dy*j +ddy;
}
