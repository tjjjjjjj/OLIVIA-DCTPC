#include "MaxCamWIMP.hh"

#include "TRandom.h"

#include "TMath.h"

#include <iostream>
#include <fstream>
#include <math.h>
using std::cout;
using std::endl;
using std::cerr;


//NamespaceImp(MaxCamWIMP)


double
MaxCamWIMP::r(double MD, double MT) {
    return 4*MD*MT/((MD+MT)*(MD+MT));
}

double
MaxCamWIMP::muSquared(double MD, double MT) {
    return 0.25*MD*MT*r(MD,MT);
}


double
MaxCamWIMP::R0(double MD, double MT, double sigma, double rhoD, double v0) {
    // Event rate per kg*day (tru), neglecting rotation around Sun and using
    // infinite escape velocity:
    //
    // MD=mass of dark matter in GeV
    // MT=mass of target nucleus in GeV
    // sigma=cross section in pb
    // rhoD=dark matter density in GeV/cm^3
    // v0=Earth velocity through dark matter halo in km/s
    //    (from galaxy rotation)
    //
    
    return 503./(MD*MT) * sigma * rhoD/0.4 * v0/230;
}



double
MaxCamWIMP::dRdE(double ER, double MD, double MT, double sigma, double rhoD, double v0, double c1, double c2) {
    // Energy distribution of recoil tracks per keV*kg*day (dru)
    // assuming infinite excape velocity and average dark matter velocity
    //

    ER*=1e-6;
    double beta=v0/3e5;
    double E0=0.5*MD*beta*beta;
    double r0=R0(MD,MT,sigma,rhoD,v0);
    double _r=r(MD,MT);
    //cout << -c2/(E0*_r) << endl;
    return c1*r0*exp(-c2*ER/(E0*_r))/(E0*_r);
}


double
MaxCamWIMP::dRdE_ERange(double E1, double E2, double MD, double MT, double sigma, double rhoD, double v0, double c1, double c2) {
    // Total rate for recoil tracks per kg*day in specific energy range (iru)
    // using the same assumptions as for dRdE.
    //

    E1*=1e-6;
    E2*=1e-6;
    double beta=v0/3e5;
    double E0=0.5*MD*beta*beta*1e6; // MD in GeV  -> E0  in keV
    double r0=R0(MD,MT,sigma,rhoD,v0);
    double _r=r(MD,MT);
    
    return r0*c1/c2*( exp(-c2*E1/(E0*_r)) - exp(-c2*E2/(E0*_r)) );
}






double
MaxCamWIMP::Maxwell(TVector3 v, TVector3 vE, double v0) {
    // Assume Maxwell-Boltzman distribution for WIMPS.
    // All velocities are in units km/s
    //

    double x=(v-vE).Mag2()/(v0*v0);
    return exp(-x);
}


void
MaxCamWIMP::getCygnusDirection(TVector3 &vec, float timeOfDay) {
    // this is approximation of Cygnus direction that assumes
    // it's on circle with diameter 90deg, centered at (0, 0) deg

    //double lat0 = 0;
    //double lon0 = 0;
    double rphi = 0.785; // 45deg
    double cosrphi=cos(rphi);
    double sinrphi=sqrt(1-cosrphi*cosrphi);
    
    double phi = timeOfDay*TMath::TwoPi();    

    vec = TVector3( cosrphi, sinrphi*cos(phi), sinrphi*sin(phi) );
}


void
MaxCamWIMP::generateMaxwellWIMP(TVector3 &vec, TVector3 vE, double v0) {
    // Generate the velocity for a WIMP based on Maxwell distribution.
    // No yearly modulation.
    
    double vx,vy,vz;
    double vEsc=800*2; // limit for velocity
    
    double rnd, fmax;
    while (1) {
        vx=(gRandom->Uniform()-0.5)*vEsc;
        vy=(gRandom->Uniform()-0.5)*vEsc;
        vz=(gRandom->Uniform()-0.5)*vEsc;
        vec=TVector3(vx,vy,vz);
        
        rnd=gRandom->Rndm();
        fmax=Maxwell( vec, vE, v0);
        //cout << rnd << "  " << fmax<<endl;
        if (fmax>rnd) break;
    }
}


double
MaxCamWIMP::nuclearFormFactor(double q)  { q++; return 1; }


double
MaxCamWIMP::CwSquared(int Z, int A) {
    // WIMP-nucleon spin factor.
    // Use  u,d: 2/3, -1/3 and 0.83, -0.43
    
    int dN = Z   + (A-Z)*2;
    int uN = Z*2 + (A-Z);
    double cw = uN*4./9.*0.83 - dN/9.*0.43;
    
    return cw*cw;
}
