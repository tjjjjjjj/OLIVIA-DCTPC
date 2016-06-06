#include "Fit/Fitter.h" 
#include "Fit/UnBinData.h" 
#define HAVE_MATHMORE
#include "DmtpcMath.hh" 
#include "TCanvas.h" 
#include "TF1.h" 
#include "MaxCamSRIM.hh" 
#include <iostream>
#include "TApplication.h" 


int main (int nargs, char ** args)
{
  TApplication app("app",0,0); 

  MaxCamSRIM srim("SRIM_He_in_CF4_100Torr"); 
  srim.setPressure(75); 


  TGraph * evr = srim.getStoppingVsRange(false); 
  TGraph * evr2 = (TGraph*) srim.getStoppingVsRange(false)->Clone("clone"); ; 

  //use only first 1000 points if more than 1000 

  int n = evr->GetN() > 1000 ? 1000 : evr->GetN(); 
  DmtpcMath::Rational rat(5,5); 
  rat.InitFitParameters(n, evr->GetX(), evr->GetY(),true,5000); 

  TCanvas c; 
  c.Divide(2,1); 
  c.cd(1); 
  TF1 f("func",rat,0,evr->GetX()[n-1], rat.NPar());
  f.SetParameters(rat.Parameters()); 
  evr->Draw("alp"); 
  f.Draw("same"); 


  c.cd(2); 
  // c.cd(3); 
  evr2->Fit(&f,"VM",0, evr->GetX()[n-1]); 
  evr2->Draw("alp"); 
  c.Show(); 



  app.Run(); 

}


