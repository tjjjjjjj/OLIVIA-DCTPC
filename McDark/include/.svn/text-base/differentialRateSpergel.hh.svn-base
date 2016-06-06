
#include "TROOT.h"
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TF12.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "bookkeeping.hh"
#include "math.h"

#include <iostream>
using std::cout;
using std::endl;

Double_t dRdERbetweenVearthandInfinity2DOverR0(Double_t *x, Double_t *par) {
  // [0] == M_dark
  // [1] == A
  // [2] == v_0
  // [3] == v_E_0
  // [4] == delta_v_E_0
  // [5] == y

  // constants
  Double_t M_p=0.938271996;
  Double_t pi=TMath::Pi();
  
  // variables
  Double_t E_R=x[0];
  Double_t cosPsi=x[1];
  
  // user input parameters
  Double_t M_dark=par[0];
  Double_t A=par[1];
  Double_t v_0=par[2];
  Double_t v_E_0=par[3];
  Double_t delta_v_E_0=par[4];
  Double_t y=par[5];

  Double_t temp = 1/v_0;
  temp /= v_0;
  temp /= M_dark;
  temp /= (1./2.);
  // have to use temp, otherwise E_0 is too small

  // calculated quantities
  Double_t M_target = M_p*A;
  Double_t r = 4*M_dark*M_target/(TMath::Power((M_dark+M_target),2));
  Double_t v_min = TMath::Sqrt( E_R*temp/(r) )*v_0;
  Double_t v_earth = v_E_0 + delta_v_E_0*TMath::Sin(2*pi*y);
  
  Double_t exponent=( v_earth*cosPsi - v_min );
  exponent *= ( v_earth*cosPsi - v_min );
  exponent /= v_0;
  exponent /= v_0;
  
  Double_t value = TMath::Exp( -exponent );
  value *= ( 1. )/( 2.);
  value *= temp;
  value *= 1/r;
  
  return value;
}

Double_t dRdERbetweenVearthandInfinity2D(Double_t *x, Double_t *par) {
  // [0] == M_dark
  // [1] == A
  // [2] == v_0
  // [3] == v_E_0
  // [4] == delta_v_E_0
  // [5] == y
  // [6] == sigma_0
  // [7] == rho_D
   
  // user input parameters
  Double_t M_dark=par[0];
  Double_t A=par[1];
  Double_t v_0=par[2];
  Double_t sigma_0=par[6];
  Double_t rho_D=par[7];
   
  Double_t value = ConvertXsectToRate( ( sigma_0*dRdERbetweenVearthandInfinity2DOverR0( x , par ) ) , v_0, A, M_dark, rho_D );
   
  return value;
}

Double_t dRdERbetweenVearthandInfinity2DERProjection(Double_t *x, Double_t *par) {
  // [0] == M_dark
  // [1] == A
  // [2] == v_0
  // [3] == v_E_0
  // [4] == delta_v_E_0
  // [5] == y
  // [6] == sigma_0
  // [7] == rho_D
  // [8] == lowerERLimit
  // [9] == upperERLimit

  // variables
  Double_t E_R=x[0];
  //Double_t cosPsi=x[1]; 

  // user input parameters
  Double_t M_dark=par[0];
  Double_t A=par[1];
  Double_t v_0=par[2];
  Double_t v_E_0=par[3];
  Double_t delta_v_E_0=par[4];
  Double_t y=par[5];
  Double_t sigma_0=par[6];
  Double_t rho_D=par[7];
  Double_t lowerERLimit=par[8];
  Double_t upperERLimit=par[9];

  TF2 *tempFcn = new TF2("tempFcn",dRdERbetweenVearthandInfinity2D,lowerERLimit,upperERLimit,-1.0,1.0,10);
  // [0] == M_dark
  // [1] == A
  // [2] == v_0
  // [3] == v_E_0
  // [4] == delta_v_E_0
  // [5] == y
  // [6] == sigma_0
  // [7] == rho_D
  // [8] == lowerERLimit
  // [9] == upperERLimit
  
  tempFcn->SetParameters(M_dark,A,v_0,v_E_0,delta_v_E_0,y,sigma_0,rho_D);
  
  // so we've got the 2D function.  Now, for a given E_R(=x), project out cosPsi(=y)
  TF12 *cosPsiGivenER = new TF12("f12",tempFcn,E_R,"y");
  
  // now integrate over all values of cosPsi
  Double_t value = cosPsiGivenER->Integral( -1.0 , 1.0 );
  
  return value;
}

void differentialRateSpergel(void){
  
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // various variables to configure plotting
  const Double_t lowerERLimit=0./1.e6;
  const Double_t upperERLimit=25./1.e6;
  const Double_t upperCosPsiLimit=1.;
  const Double_t lowerCosPsiLimit=-1.;
  //const Double_t lowerXsectLimit=0;
  //const Double_t upperXsectLimit=1.;
  
  TF2 *differentialRateSpergel = new TF2("dRdERbetweenVearthandInfinity2DERProjection",
                                         dRdERbetweenVearthandInfinity2D,lowerERLimit,upperERLimit,lowerCosPsiLimit,upperCosPsiLimit,10);
  // [0] == M_dark
  // [1] == A
  // [2] == v_0
  // [3] == v_E_0
  // [4] == delta_v_E_0
  // [5] == y
  // [6] == sigma_0
  // [7] == rho_D

  // for the Earth limits, it's nice to have an index of months
  /*const Double_t daysInaYear=365.242;
  const Double_t yJanuary=-61./daysInaYear;
  const Double_t yFebruary=-30./daysInaYear;
  const Double_t yMarch=-1./daysInaYear;
  const Double_t yApril=30./daysInaYear;
  const Double_t yMay=60./daysInaYear;
  const Double_t yJune=91./daysInaYear;
  const Double_t yJuly=121./daysInaYear;
  const Double_t yAugust=152./daysInaYear;
  const Double_t ySeptember=183./daysInaYear;
  const Double_t yOctober=213./daysInaYear;
  const Double_t yNovember=244./daysInaYear;
  const Double_t yDecember=274./daysInaYear;
  */
  const Double_t c = 3.e8;
  Double_t M_dark=60.;
  Double_t A=132.;
  Double_t v_0=230.e3/c;
  Double_t rho_D=0.4;  // GeV per cm^3 c^2
  Double_t sigma_0=1e-44; // cm^2

  Double_t v_E_0=244.e3/c;
  Double_t delta_v_E_0=15.e3/c;
  Double_t y=0;
  
  differentialRateSpergel->SetParameters(M_dark,A,v_0,v_E_0,delta_v_E_0,y,sigma_0,rho_D);

  // debug scan of xsect values
  Double_t number_of_points_to_scan = 20.;
  // the +2 is to keep from running into the edges
  Double_t step_size = ( fabs( upperERLimit - lowerERLimit ) )/( number_of_points_to_scan + 2. ); 
  cout << step_size << endl;
  
  for(Int_t i = 0; i<20; i++){
    cout << "differentialRateSpergel->Eval(" << (i+1)*step_size << " , 0.5 )=" << differentialRateSpergel->Eval( (i+1)*step_size , 0.5 ) << endl;
  }

  TCanvas* c1=new TCanvas("c1","myCanvas");
  gStyle->SetCanvasColor(0);
  c1->SetFrameFillColor(0);
  //  c1->SetLogy();
  //  c1->SetLogx();
  //  differentialRate1->SetMinimum(lowerXsectLimit);
  //  differentialRate1->SetMaximum(upperXsectLimit);
  differentialRateSpergel->SetTitle("");
  differentialRateSpergel->GetXaxis()->SetTitle("E_{R} [GeV]");
  differentialRateSpergel->GetYaxis()->SetTitle("Cos(#psi)");
  differentialRateSpergel->GetXaxis()->SetTitleOffset(2.);
  differentialRateSpergel->GetYaxis()->SetTitleOffset(2.);
  differentialRateSpergel->GetXaxis()->CenterTitle();
  differentialRateSpergel->GetYaxis()->CenterTitle();

  //  leg = new TLegend(0.48,0.64,0.83,0.88);  //coordinates are fractions
  //  leg->AddEntry(differentialRate1,"0<v<#infty","l");  // "l" means line
  //  leg->AddEntry(differentialRate2,"0<v<v_{esc}","l");  // "l" means line
  //  leg->AddEntry(differentialRate3,"v_{Earth}<v<#infty","l");  // "l" means line
  //  leg->AddEntry(differentialRateSpergel,"v_{Earth}<v<v_{esc}","l");  // "l" means line
  //  leg->AddEntry(differentialRateLS,"v_{Earth}<v<#infty (LS param)","l");
  //  leg->SetBorderSize(0);                                         //of pad dimensions
  //  leg->SetTextSize(0.05);
  //  leg->SetFillColor(0);
  
  differentialRateSpergel->Draw("LEGO");
  //  leg->Draw();

  //  c1->SaveAs("differentialRateComparison.eps");
} 

