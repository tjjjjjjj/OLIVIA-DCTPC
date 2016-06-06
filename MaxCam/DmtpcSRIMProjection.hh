#ifndef DMTPC_SRIM_PROJECTION_HH
#define DMTPC_SRIM_PROJECTION_HH

#include "TGraph.h"
#include "TSpline.h"
#include "GausConvSpline.hh"
#include "TF1.h" 

class MaxCamSRIM; 

class DmtpcSRIMProjection 
{

  public: 

    enum PARTICLE { ALPHA, F, C}; 
    DmtpcSRIMProjection(PARTICLE type, double blur_amt, double pressure, double a,  double b); 
    DmtpcSRIMProjection(MaxCamSRIM *srim, double blur_amt, double pressure, double a,  double b); 
    virtual ~DmtpcSRIMProjection(); 
    double dE(double x, double nsig=5) const; 
    double Eval(double *x, double *p) {return dE(*x,*p); }
    double E(double x0, double x1, double nsig = 5); 
    void setGain(double adu_per_kev) { gain = adu_per_kev; }
    void setSinTheta(double sin_theta) { this->sin_theta = sin_theta; }  
    void setOffset(double offset) { this->offset=offset; }
    void setLengthCal(double mm_per_px) { lengthcal = mm_per_px;}  
    TF1 * asFunction(const char * name, double nsig_bounds = 5); 
    const GausConvSpline * spline() const { return stoppingR; } 
    const TSpline3 * rawSpline() const { return srim_stoppingRSpline; } 
//    double totalEatPosition(double x); 
//
    static const char * getFileName(PARTICLE type); 


 private: 

    GausConvSpline * stoppingR; 
    TSpline3 * srim_stoppingRSpline; 
    double a; 
    double b; 
    double blur; 
    double gain; 
    double offset; 
    double lengthcal; 
    double sin_theta; 
    TF1 * func; 


 ClassDef(DmtpcSRIMProjection,0); 
}; 


#endif 
