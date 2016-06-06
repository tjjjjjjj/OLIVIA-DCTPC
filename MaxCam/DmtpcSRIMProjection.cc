#include "DmtpcSRIMProjection.hh"
#include "MaxCamSRIM.hh"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/Integrator.h"
#include <iostream>

ClassImp(DmtpcSRIMProjection); 



TF1 * DmtpcSRIMProjection::asFunction(const char * name, double nsig_bounds) 
{
  TF1 * func = new TF1(name,this,&DmtpcSRIMProjection::Eval, a - blur * nsig_bounds, b + blur * nsig_bounds, 
                                                      1, "DmtpcSRIMProjection","Eval"); 
  func->SetParameter(0, 5); 
  return func; 
}


const char * DmtpcSRIMProjection::getFileName(PARTICLE type)
{
  switch(type)
  {
    case ALPHA:               
      return "SRIM_He_in_CF4_100Torr"; 
    case C: 
      return "SRIM_C_in_CF4_100Torr"; 
    case F: 
      return "SRIM_F_in_CF4_100Torr"; 

    default: 
      std::cerr << "Invalid Type." << std::endl; 
  }

  return 0; 

}


double DmtpcSRIMProjection::dE(double x, double nsig) const
{

  //convert to physical units along track 

  x-= offset; 
  x*= lengthcal;  //pixel to mm 
  x/= sin_theta;  //convert to track direction
  double val = stoppingR->Eval(&x,&nsig); 
  //convert to adu 
  return val * gain; 
}

//possible to do this analytically I think, but not trying... 
double DmtpcSRIMProjection::E(double x0, double x1, double nsig)
{

  ROOT::Math::WrappedTF1 wf1(*func);
  wf1.SetParameters(&nsig); 

  ROOT::Math::Integrator ig; 
      
  ig.SetFunction(wf1, false);
  ig.SetRelTolerance(0.001);

  return ig.Integral(x0,x1); 

}


DmtpcSRIMProjection::DmtpcSRIMProjection(MaxCamSRIM * sr, double blur_amt, double pressure, double a, double b) 
{

  this->a = a; 
  this->b = b; 
  this->blur = blur_amt; 

  this->gain = 1; 
  this->lengthcal = 1; 
  this->sin_theta = 1; 


  sr->setPressure(pressure); 
  TGraph * srim_stopping_r = sr->getStoppingVsRange(false); 
  srim_stoppingRSpline = new TSpline3("org_spline", srim_stopping_r); 
  stoppingR = new GausConvSpline(blur_amt, srim_stoppingRSpline, a, b); 

  func = asFunction("srim_proj_func"); 
//  delete srim_stopping_r; 
}



DmtpcSRIMProjection::DmtpcSRIMProjection(PARTICLE type, double blur_amt, double pressure, double a, double b) 
{
  const char * srim_str = getFileName(type); 

  this->a = a; 
  this->b = b; 
  this->blur = blur_amt; 

  this->gain = 1; 
  this->lengthcal = 1; 
  this->sin_theta = 1; 


  MaxCamSRIM sr(srim_str); 
  sr.setPressure(pressure); 
  TGraph * srim_stopping_r = sr.getStoppingVsRange(false); 
  srim_stoppingRSpline = new TSpline3("org_spline", srim_stopping_r); 
  stoppingR = new GausConvSpline(blur_amt, srim_stoppingRSpline, a, b); 

  func = asFunction("srim_proj_func"); 
//  delete srim_stopping_r; 
}


DmtpcSRIMProjection::~DmtpcSRIMProjection()
{
  delete srim_stoppingRSpline; 
  delete stoppingR; 
  delete func; 
}
