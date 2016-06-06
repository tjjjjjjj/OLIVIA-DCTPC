#ifndef MAXCAM_TWOBODY_KINEMATICS_HH
#define MAXCAM_TWOBODY_KINEMATICS_HH

#include <assert.h>

#include "TROOT.h"

namespace MaxCamTwoBodyKinematics {

    double calcCosRecoilFromRecoilEnergy(double P, double p0, double M, double m);
    double calcRecoilEnergyFromCosRecoil(double cosTheta, double p0, double M, double m);
    
    double calcRecoilEnergyFromCosScatterCMS(double cosTheta, double p0, double M, double m);
    double calcCosScatterCMSFromRecoilEnergy(double P, double p0, double M, double m);
    
    double calcCosScatterFromRecoilEnergy(double E_M, double E0, double M, double m);
    
    double Legendre(double x, double *par, int n);

};

#endif

