#include "MaxCamTwoBodyKinematics.hh"

#include <iostream>
#include <fstream>
#include <math.h>
using std::cout;
using std::endl;
using std::cerr;


//NamespaceImp(MaxCamTwoBodyKinematics)



double
MaxCamTwoBodyKinematics::calcCosRecoilFromRecoilEnergy(double Er, double E0, double M, double m) {
    // Elastic, non-relativistic two body scattering.
    // Computes cosine of recoil angle in lab frame for given recoil 
    // energy (Er), mass (M) of recoil particle, initial energy of
    // projectile (E0) and mass of projectile (m):
    // 
    //  cos(Theta) = 1/2 * sqrt(Er/E0) * (1+M/m) / sqrt(M/m)
    //
    //                     (E,M) 'recoil'
    //    'projectile'               
    //    (E0,m) ---> (0,M) 'target'  
    //                   
    //                     (p,m) 'scatter'
    //
    // (inverse of calcRecoilEnergyFromCosRecoil)
    
    double massRatio = M/m; // target over projectile
    return sqrt(Er/(massRatio*E0))*(1+massRatio)*0.5;        
}

double
MaxCamTwoBodyKinematics::calcRecoilEnergyFromCosRecoil(double cosTheta, double E0, double M, double m) {
        // Elastic, non-relativistic two body scattering.
        // Computes recoil energy in lab frame for given recoil
        // angle (cosTheta), mass (M) of recoil particle, initial energy of projectile (E0)
        // and mass of projectile (m):
        // 
        //
        //  Er = 4 * cosTheta^2 * M/m * E0 / (1+M/m)^2
        //
        // (inverse of calcCosRecoilFromRecoilEnergy)

        double massRatio = M/m; // target over projectile
        return 4*cosTheta*cosTheta*massRatio*E0/((1+massRatio)*(1+massRatio));
}


double 
MaxCamTwoBodyKinematics::calcRecoilEnergyFromCosScatterCMS(double cosTheta, double E0, double M, double m) {
        // Elastic, non-relativistic two body scattering.
        // Computes recoil energy in lab frame for given projectile scattering angle
        // in CMS (cosTheta), mass (M) of recoil particle, initial energy (E0)
        // and mass (m) of projectile:
        // 
        //  P = (2*E0*M/m) / (1+M/m)^2 * (1-cosTheta)
        //
        // (inverse of calcCosScatterFromRecoilEnergy)

  double massRatio = M/m; // target over projectile
  return (E0*2*massRatio) / ((1+massRatio)*(1+massRatio)) * (1-cosTheta);
}


double
MaxCamTwoBodyKinematics::calcCosScatterCMSFromRecoilEnergy(double Er, double E0, double M, double m) {
        // Elastic, non-relativistic two body scattering.
        // Computes scattering angle of projectile in CMS frame for given recoil energy 
        // in LAB frame (Er), mass (M) of recoil particle, initial energy (E0)
        // and mass (m) of projectile:
        // 
        //  cosScatterCMS = 1 - Er*(1+M/m)^2/(2*p0*M/m)
        //
        // (inverse of calcRecoilEnergyFromScatterAngle)

  double massRatio = M/m; // target over projectile
  return 1 - Er*(1+massRatio)*(1+massRatio)/(2*massRatio*E0);
}



double
MaxCamTwoBodyKinematics::calcCosScatterFromRecoilEnergy(double E_M, double E0, double M, double m) {
    // Elastic, non-relativistic two body scattering.
    // Computes scattering angle of projectile (m) in lab frame for given recoil energy
    // of mass (M), initial energy (E0):
    // 
    //      sin_m * p_m = sin_M * p_M
    //

    double cos_M = calcCosRecoilFromRecoilEnergy( E_M, E0, M, m);
    double sin_M = sqrt( 1 - cos_M*cos_M );

    double E_m = E0-E_M; 
    double p_m = sqrt(2*m*E_m);
    double p_M = sqrt(2*M*E_M);
    double sin_m = p_M * sin_M / p_m;
    
    return sqrt( 1 - sin_m*sin_m );
}


double 
MaxCamTwoBodyKinematics::Legendre(double x, double *par, int n) {
        // Legendre polynomials for angular distributions of 
        // scattering products.

        assert(par);
        if (n<1) return 1;
        double *p=new double[n];
        p[0]=1; 
        p[1]=x; 
        for (int i=1; i<n-1; i++) {
            p[i+1] = ( (2*i+1)*x*p[i] - i*p[i-1] )/( i+1 );
        }
        double sum = 0;
        for (int i=0; i<n; i++) sum += 0.5*(2*i+1) * p[i] * par[i];
        delete p;
        return sum;
}
