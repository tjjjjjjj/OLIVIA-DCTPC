// code that spits out the functional form of the differential rate
// as a function of nuclear kinetic recoil energy
#include "differentialRateSpergel.hh"


// code that loads root stuff
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h" 
#include "TH2F.h"

// nEvt is the number of E_R's to throw
// where E_R is the nuclear kinetic recoil energy
// if debug is >0, lots of printouts.  Else, quiet
// freq is the frequency of printouts
void particleGun2D(Int_t nEvt=1, Int_t debug=0, Int_t freq=1.e5){
  
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // various variables to configure plotting
  // the differential rate will be computed between these two limits of nuclear kinetic recoil energy, in GeV.  
  const Double_t lowerERLimit=0./1.e6; // to put into GeV
  const Double_t upperERLimit=25./1.e6; // to put into GeV
  const Double_t lowerCosPsiLimit=-1.0;
  const Double_t upperCosPsiLimit=1.0;

  // histogram to hold all of the E_R's we're going to pop
  TH2F *hdNdE;
  hdNdE = new TH2F("hdNdE", "", 75,lowerERLimit,upperERLimit,75,lowerCosPsiLimit,upperCosPsiLimit);

  TF2 *differentialRateSpergel = new TF2("dRdERbetweenVearthandInfinity2DERProjection",dRdERbetweenVearthandInfinity2D,lowerERLimit,upperERLimit,lowerCosPsiLimit,upperCosPsiLimit,10);
  // [0] == M_dark
  // [1] == A
  // [2] == v_0
  // [3] == v_E_0
  // [4] == delta_v_E_0
  // [5] == y
  // [6] == sigma_0
  // [7] == rho_D

  // the differential rate has a dependence on where the Earth is around the sun; the following are some useful values of y
  const Double_t daysInaYear=365.242;
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
  Double_t v_esc=600.e3/c;

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

  // debug scan of differential rate values, to make sure they seem reasonable.
  if( debug>0 ){
    Double_t number_of_points_to_scan = 20.;
    // the +2 is to keep from running into the edges
    Double_t step_size = ( fabs( upperERLimit - lowerERLimit ) )/( number_of_points_to_scan + 2. ); 
    cout << step_size << endl;
    
//    for(Int_t i = 0; i<20; i++){
//      cout << "differentialRate4->Eval(" << (i+1)*step_size << ",0.5)=" << differentialRate4->Eval( (i+1)*step_size , 0.5 ) << endl;
//    }
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
  differentialRateSpergel->GetZaxis()->SetTitle("d^{2}R/dE_{R}d(Cos(#psi)) [GeV^{-1}kg^{-1}day^{-1}]");

  differentialRateSpergel->GetXaxis()->SetTitleOffset(2.);
  differentialRateSpergel->GetYaxis()->SetTitleOffset(2.);
  differentialRateSpergel->GetZaxis()->SetTitleOffset(1.3);
  differentialRateSpergel->GetXaxis()->CenterTitle();
  differentialRateSpergel->GetYaxis()->CenterTitle();
 
  differentialRateSpergel->Draw("LEGO");
  c1->SaveAs("differentialRateSpergel_ParticleGun2D.eps");

  for( Int_t i=1; i<(nEvt+1); i++ ){
    
    // note that you get the energies from this back in GeV, although
    // keV's are the more natural energy scale in this problem
    // E_R is what you want to use in the particle gun; pull it randomly from the differential rate
    Double_t E_R=-999.;
    Double_t cosPsi=-999.;
    differentialRateSpergel->GetRandom2(E_R,cosPsi);
    hdNdE->Fill( E_R , cosPsi );
    if(debug>0 || (i+1)%freq==0){ 
      cout << "**DEBUG event " << (i+1) 
	   << " E_R=" << E_R*1.e6 << " keV" 
	   << " Cos(Psi)=" << cosPsi << endl; 
    }  
    
    
  }

  TCanvas* c2=new TCanvas("c2","myCanvas");
  hdNdE->SetXTitle("E_{R} [GeV]");
  hdNdE->SetYTitle("Cos(#psi)");
  hdNdE->SetZTitle("d^{2}N/dE_{R}d(Cos(#psi))");
  hdNdE->GetXaxis()->SetTitleOffset(2.);
  hdNdE->GetYaxis()->SetTitleOffset(2.);
  hdNdE->GetZaxis()->SetTitleOffset(1.1);
  hdNdE->Draw("LEGO");
  c2->SaveAs("dNdE_ParticleGun2D.eps");
}

