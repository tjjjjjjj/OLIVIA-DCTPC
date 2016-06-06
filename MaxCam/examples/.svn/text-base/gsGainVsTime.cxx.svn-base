/************************************************************
//
//  A script for analyzing data.
//  Run the script on any machine from cint:
//
//  .L gainVsTime.cxx
//  all()
//
*************************************************************/

#include <ctime>
#include <cmath>
#include <sstream>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <map>
// 
// 
// //#include "TFile.h"
// //#include "TChain.h"
// //#include "TLeaf.h"
// 
// 
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
// 
#include "../MaxCamRead.hh"
#include "../MaxCamChannel.hh"
// 
// using namespace std;
// 

MaxCamRead *ana;
int i, n, irej;       //i=envtNumber; n=NumberOfEvents, irej=#rejected events

TF1 *fun;
float x[9];           // y projection of picture: time,y,yerr,width,chi2,ndf,wirehv,meshhv,pressure
long time0;
TGraphErrors *gr=0;
float sqrt2pi=sqrt(2*TMath::Pi());
float ccdGain=1.7;   // every PE is 1.7 counts from Flat Field calibration - TBC 
TTree *nt;


struct wireData {
  float y, yerr, width, chi2;
  int ndof;
} wire0, wire1, wire2, wire3, wire4, wire5, *wire;

void init(int runN) {

  //  ana = new MaxCamRead("/data/ccdrun_00045.root");
  //  ana = new MaxCamRead("../../Runs/ccdrun_00052.root"); // wire1=bin4-13; 2=23-27; 3=40-45; 

  char filename[256];
  snprintf(filename,sizeof(filename),"%s%d%s","../../Runs/ccdrun_000",runN,".root"); 
  cout << " Filename = " << filename << endl; 
  if(runN !=0) { 
    ana = new MaxCamRead(filename);
    ana->readBiasFrame(); // reads bias frame and automatically subtracts it from all events 
                        // to stop this: deleteBiasFrame()

  } else { 
    // all all good and compatible runs to increase stat ... 
    ana = new MaxCamRead("../../Runs/ccdrun_00048.root");
    ana->readBiasFrame(); // reads bias frame and automatically subtracts it from all events 
                        // to stop this: deleteBiasFrame()  
    ana->addRun("../../Runs/ccdrun_00049.root"); 
    ana->addRun("../../Runs/ccdrun_00050.root"); 
    ana->addRun("../../Runs/ccdrun_00051.root"); 
    ana->addRun("../../Runs/ccdrun_00052.root"); 
    ana->addRun("../../Runs/ccdrun_00053.root"); 
    //    ana->addRun("../../Runs/ccdrun_00055.root"); // not triggered by eye 
    //  ana->addRun("../../Runs/ccdrun_00056.root");   // mesh scan 
  }

  i=0; 
  irej=0;

  fun = new TF1("fun","[0]+[1]*exp(-0.5*(x-[2])**2/[3]**2)"); 

  n=ana->tree()->GetEntries(); 
  cout << "TOTAL = " << n << endl;

  ana->getEvent(0);
  time0=ana->timeStamp()->Get();

  nt= new TTree("res","");
  //  nt->Branch("fit", &x[0], "time/F:y:yerr:width:chi2:ndf:wirehv:meshhv:pressure");
  nt->Branch("conf",  &x[0],  "time/F:wirehv:meshhv:pressure");
  nt->Branch("wire0", &wire0, "y/F:yerr:width:chi2:ndof/I");
  nt->Branch("wire1", &wire1, "y/F:yerr:width:chi2:ndof/I");
  nt->Branch("wire2", &wire2, "y/F:yerr:width:chi2:ndof/I");
  nt->Branch("wire3", &wire3, "y/F:yerr:width:chi2:ndof/I");
  nt->Branch("wire4", &wire4, "y/F:yerr:width:chi2:ndof/I");
  nt->Branch("wire5", &wire5, "y/F:yerr:width:chi2:ndof/I");

}

void eventConfig() {
  x[0] = ana->timeStamp()->Get()-time0;
  x[1] = ana->wire()->currentValue;
  x[2] = ana->mesh()->currentValue;
  x[3] = ana->pressure()->currentValue;
}

void wireYield(int iwire,int runNumber ) { 
  
  // default: case 0
  int wireMinPixel=4; 
  int wireMaxPixel=10; // to be fixed!   
  wire=&wire0;
  
  
  if( runNumber == 52 || runNumber == 53 ||runNumber == 55 ||runNumber == 57 || runNumber == 0  ) { 
    switch(iwire) {
    case 1: wireMinPixel=23; wireMaxPixel=27; wire=&wire1; break;
    case 2: wireMinPixel=40; wireMaxPixel=44; wire=&wire2; break;
    case 3: wireMinPixel=58; wireMaxPixel=62; wire=&wire3; break;
    case 4: wireMinPixel=73; wireMaxPixel=76; wire=&wire4; break;
    case 5: wireMinPixel=93; wireMaxPixel=96; wire=&wire5; break;
    }
    
  } else {
    cout << " ======== Problem: find position of wires! ===== " << endl; 
  }
  
//   cout << "Run number="<< runNumber << " select wireNumber " << iwire
//        << " in  bins " <<wireMinPixel <<"-"<< wireMaxPixel << endl; 
  
  
  
  // find intensity
  TH1D *hy = ana->ccdImage()->ProjectionY("_py",wireMinPixel,wireMaxPixel,"e"); // use only one wire
  double brho = hy->GetNbinsX()/(hy->GetXaxis()->GetXmax()-hy->GetXaxis()->GetXmin());
  
  fun->SetParameters( hy->GetMinimum(),
		      hy->GetMaximum()-hy->GetMinimum(), 
		      hy->GetBinCenter(hy->GetMaximumBin()), // educated guesses
		      5);


  hy->Fit("fun","Q0");
  //c1->Update(); gSystem->Sleep(1000); 


  float I = hy->GetFunction("fun")->GetParameter(1)
    * hy->GetFunction("fun")->GetParameter(3)
    * sqrt2pi*brho
    / ccdGain;

  float IErr = I* sqrt(
    pow(hy->GetFunction("fun")->GetParError(1)/hy->GetFunction("fun")->GetParameter(1),2) +
    pow(hy->GetFunction("fun")->GetParError(3)/hy->GetFunction("fun")->GetParameter(3),2) );

  // time stamp
  wire->y     = I; 
  wire->yerr  = IErr;
  wire->width = hy->GetFunction("fun")->GetParameter(3);
  wire->chi2  = hy->GetFunction("fun")->GetChisquare();
  wire->ndof  = hy->GetFunction("fun")->GetNDF();

  //  cout << wire->y << endl;

}


// bool wireYield(int wireMinPixel=1, int wireMaxPixel=96) { // arguments: bins of interest for a wire in X 

//   //  int ii=i-irej;
//   ana->getEvent(i++);
  
//   // find intensity
//   TH1D *hy = ana->ccdImage()->ProjectionY("_py",wireMinPixel,wireMaxPixel,"e"); // use only one wire
//   double brho = hy->GetNbinsX()/(hy->GetXaxis()->GetXmax()-hy->GetXaxis()->GetXmin());// bin density 
  
//   fun->SetParameters( hy->GetMinimum(),// P0=flat background 
// 		      hy->GetMaximum()-hy->GetMinimum(), // P1=height of gaussian 
// 		      hy->GetBinCenter(hy->GetMaximumBin()), // P2=mean : educated guesses
// 		      5.);// P3= width 


//   hy->Fit("fun","Q0"); // Q=quiet fit; 0=no graphic output 
//   //  c1->Update(); gSystem->Sleep(1000); // if you want to look at the events 


//   // Intensity=integral of gaussian / CCD gain 
//   float I = hy->GetFunction("fun")->GetParameter(1) 
//     * hy->GetFunction("fun")->GetParameter(3)
//     * sqrt2pi*brho
//     / ccdGain;

//   float IErr = I* sqrt(
//     pow(hy->GetFunction("fun")->GetParError(1)/hy->GetFunction("fun")->GetParameter(1),2) +
//     pow(hy->GetFunction("fun")->GetParError(3)/hy->GetFunction("fun")->GetParameter(3),2) );

//   // time stamp
//   x[0] = ana->timeStamp()->Get()-time0;
//   x[1] = I; 
//   x[2] = IErr;
//   x[3] = hy->GetFunction("fun")->GetParameter(3);
//   x[4] = hy->GetFunction("fun")->GetChisquare();
//   x[5] = hy->GetFunction("fun")->GetNDF();
//   x[6] = ana->wire()->currentValue;
//   x[7] = ana->mesh()->currentValue;
//   x[8] = ana->pressure()->setValue;

//   nt->Fill();

//   return true;
// }


void gsGainVsTime(int  runNumber = 0) {

  init(runNumber);

   //  while (i<=n-1) wireYield(minX,maxX);
   while (i<=n-1){
     ana->getEvent(i);
     eventConfig();
     for (int iw=0;iw<6;iw++) wireYield(iw,runNumber);
     nt->Fill();
     cout <<"event="<<i<<endl;
     i++;
   }  


  
  // plot all wire yields vs time: nt->Draw("y:time");
  // plot image for event 34: ana->getEvent(34); ana->ccdImage()->Draw("colz");

//old cuts "pressure>0&&y<25000&&y>0&&wirehv>2.0&&abs(meshhv-1.5)<0.1&&abs(pressure-190)<10" 

   // one simple example. 
//  TCanvas* canvas = new TCanvas("mycan","mycan",400,400);   canvas->Divide(1,1);  canvas->cd(1);
//   nt->Draw("log(wire0.y):conf.wirehv","conf.wirehv>2.3");
//  canvas->Print("~/public_html/tmp/test.gif"); 


// Gain vs  Wire hv for all wires in the run 
//TCanvas* canvas = new TCanvas("mycan","mycan",400,600);   canvas->Divide(2,3);  
//  canvas->cd(1); nt->Draw("(wire0.y):conf.wirehv","conf.wirehv>2.3 && wire0.y>0 && wire0<10000");
//  canvas->cd(2); nt->Draw("(wire1.y):conf.wirehv","conf.wirehv>2.3 && wire1.y>0 && wire1<10000");
//  canvas->cd(3); nt->Draw("(wire2.y):conf.wirehv","conf.wirehv>2.3 && wire2.y>0 && wire2<10000");
//  canvas->cd(4); nt->Draw("(wire3.y):conf.wirehv","conf.wirehv>2.3 && wire3.y>0 && wire3<10000");
//  canvas->cd(5); nt->Draw("(wire4.y):conf.wirehv","conf.wirehv>2.3 && wire4.y>0 && wire4<10000");
//  canvas->cd(6); nt->Draw("(wire5.y):conf.wirehv","conf.wirehv>2.3 && wire5.y>0 && wire5<10000");
//  //canvas->Print("~/public_html/tmp/Gain_vs_HV_6wires.gif"); 
//  canvas->Print("~/public_html/tmp/Run52_Gain_vs_HV_6wires.gif"); 
//
//TCanvas* canvas = new TCanvas("mycan","mycan",400,400);   canvas->Divide(1,1);  
//  canvas->cd(1); nt->Draw("(wire0.y+wire1.y+wire2.y):conf.wirehv","conf.wirehv>2.3 && wire0.y>0 && wire0<10000");
//canvas->Print("~/public_html/tmp/test.gif"); 
//
//
//  TCanvas* canvas = new TCanvas("mycan","mycan",400,400);   canvas->Divide(1,1);  canvas->cd(1);
//   nt->Draw("(wire0.y)","conf.wirehv>2.3&&wire1.y<2000");                                
//   canvas->Print("~/public_html/tmp/run53_GainWire0_Wire1_lt_2000.gif");                                                 
//
//  TCanvas* canvas = new TCanvas("mycan","mycan",400,400);   canvas->Divide(1,1);  canvas->cd(1);
//   nt->Draw("(wire0.y)","conf.wirehv>2.3&&wire1.y>2000");                                       
//  canvas->Print("~/public_html/tmp/run53_GainWire0_Wire1_gt_2000.gif");                                                 
//
//

   TCanvas* canpoo = new TCanvas("canpoo","canpoo",600,600);
   canpoo->Divide(6,3);

   const int NHV = 6;
   const int NWIRE = 3;
   const double hv[NHV] = { 2.55, 2.60, 2.65, 2.70, 2.75, 2.80 };
   const double cutwire1[NHV] = { 700., 900., 1200., 1700., 2200., 2500. };
   const double cutwire0[NHV] = { 1500., 1900., 2500., 3400., 4500., 5000. };

   double gMean[NWIRE][NHV];
   double gSigma[NWIRE][NHV];

   TH1F* hGain[NHV][NWIRE];
   char hname[100];
   char vname[100];
   char cut[200];
   for (int iHV=0; iHV<NHV; iHV++){
     for (int iWIRE=0; iWIRE<NWIRE; iWIRE++) {
       sprintf(hname, "hGain_%d_%d", iHV, iWIRE);
       hGain[iHV][iWIRE] = new TH1F(hname,"Gain",100,0.,5000.);
       sprintf(vname, "wire%d.y", iWIRE);
       if (iWIRE==0) {
	 sprintf(cut, "wire1.y<%f && abs(conf.wirehv-%f)<0.020", cutwire1[iHV], hv[iHV]);
       } else {
	 sprintf(cut, "wire0.y<%f && abs(conf.wirehv-%f)<0.020", cutwire0[iHV], hv[iHV]);
       }
       nt->Project(hname, vname, cut);
       canpoo->cd(iHV+1+6*iWIRE);
       hGain[iHV][iWIRE]->Draw();
       gMean[iWIRE][iHV] = hGain[iHV][iWIRE]->GetMean();
       gSigma[iWIRE][iHV] = hGain[iHV][iWIRE]->GetMeanError();
       cout << "HV = " << hv[iHV] << ", wire " << iWIRE << ", Gain = " << gMean[iWIRE][iHV] << " +/- " << gSigma[iWIRE][iHV] << endl;
     }
   }

   //   canpoo->Print("~/public_html/tmp/canpoo.gif");

   TCanvas* cangain = new TCanvas("cangain","cangain",600,600);
   cangain->Divide(2,2);

   TH1F* gainVsHv[NWIRE];
   //   const double zero[NHV] = { 1., 1., 1., 1., 1., 1. };
   for (int iWIRE = 0; iWIRE<NWIRE; iWIRE++) {
     sprintf(hname,"gainVsHv_wire%d",iWIRE);
     gainVsHv[iWIRE] = new TH1F(hname, "Gaiv vs HV", NHV, 2.525, 2.825);
     for (int iHV=0; iHV<NHV; iHV++) {
       gainVsHv[iWIRE]->SetBinContent(iHV+1,gMean[iWIRE][iHV]);
       gainVsHv[iWIRE]->SetBinError(iHV+1,gSigma[iWIRE][iHV]);
     }
     cangain->cd(iWIRE+1);
     gainVsHv[iWIRE]->Draw();
   }

   //   cangain->Print("~/public_html/tmp/cangain.gif");

//   if(nt->abs(conf.wirehv-2.55)<0.020) { 
//   wire 0: wire1.y<700
//        1+: wire0.y<1500
//
//   abs(conf.wirehv-2.60)<0.020
//   wire 0: wire1.y<900
//        1+: wire0.y<1900
//
//   abs(conf.wirehv-2.65)<0.020
//   wire 0: wire1.y<1200
//        1+: wire0.y<2500
//
//   abs(conf.wirehv-2.7)<0.020
//   wire 0: wire1.y<1700
//        1+: wire0.y<3400
//
//
//   abs(conf.wirehv-2.75)<0.020
//   wire 0: wire1.y<2200
//        1+: wire0.y<4500
//
//   abs(conf.wirehv-2.8)<0.020
//   wire 0: wire1.y<2500
//        1+: wire0.y<5000
//
//
  
}
