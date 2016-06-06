#ifndef MAXCAM_WIMP_HH
#define MAXCAM_WIMP_HH

#include "TROOT.h"
#include "TVector3.h"


namespace MaxCamWIMP {

    double R0(double MD, double MT, double sigma=1, double rhoD=0.4, double v0=230);
    
    double dRdE(double ER,
                double MD, double MT, double sigma=1, double rhoD=0.4, double v0=230,
                double c1=0.751, double c2=0.561);

    double dRdE_ERange(double E1, double E2,
                       double MD, double MT, double sigma=1, double rhoD=0.4, double v0=230,
                       double c1=0.751, double c2=0.561);
   
    double Maxwell(TVector3 v, TVector3 vE, double v0=230);

    void generateMaxwellWIMP(TVector3 &v, TVector3 vE, double v0=230);
    
    void getCygnusDirection(TVector3 &v, float timeOfDay);

    double nuclearFormFactor(double q); 

    double CwSquared(int Z, int A);

    double r(double MD, double MT);

    double muSquared(double MD, double MT);
    
};

#endif

