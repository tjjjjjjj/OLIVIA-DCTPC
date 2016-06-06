#define DCTPCTree_cxx
#include "../input_files/DCTPCTree.h"
#include "../input_files/Calibration_constants_LittleDCTPC.h";
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
using namespace std;

void DCTPCTree::Loop()
{
TStopwatch timer;
  gROOT->SetStyle("Plain");
  gStyle->SetEndErrorSize(3);
  gStyle->SetPalette(1,0);
  gStyle->SetLineWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptStat(kTRUE);
  gStyle->SetOptFit(kTRUE);
  TH1::AddDirectory(false);
  
  TFile *outtree = new TFile("$LittleDCTPC_physics_input");
  TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
  DCTPCTree aStep(dctreepc);
  
   TH2D *hist_ccd_vs_length=new TH2D("CCD vs Length", "CCD vs Length", 100, 0, 10000, 100, 0, 10000);
  
//create event loop
  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;

  int passccd=0;
  int passwfandccd=0;

  for (int event = 0; event<nentries; event++)
    {
      aStep.GetEntry(event);
      if(event%1000000==0)
      cout<<((double)event/(double)nentries)*100.<<"%"<<endl;

      r2 = (aStep.Track_x_pix*aStep.Track_x_pix+aStep.Track_y_pix*aStep.Track_y_pix)*pow(MMPERPIXEL/10.,2);
      weight=1./pow(RADIAL_PARAM0_CALIB/(RADIAL_PARAM1_CALIB+r2),RADIAL_PARAM2_CALIB);//radial calibration 
      verticallength=((aStep.Mesh_R0time_samp+aStep.Mesh_F10time_samp)*driftspeed);
      horizontallength=aStep.Track_range_pix*MMPERPIXEL;
      tracklength=(sqrt((verticallength*verticallength)+(horizontallength*horizontallength)));
      lengthenergy=El_0 + (El_1*tracklength) + (El_2*tracklength*tracklength);
      anodecalib=ANODE_CALIB;
      meshcalib=MESH_CALIB;
      
      delta_anodeccd=((aStep.Etrig_kev*anodecalib)-(aStep.Etrack_kev*weight))/((aStep.Etrig_kev*anodecalib)+(aStep.Etrack_kev*weight));
      delta_meshccd=((aStep.Emesh_kev)-(aStep.Etrack_kev*weight))/((aStep.Emesh_kev)+(aStep.Etrack_kev*weight));
      delta_anodemesh=((aStep.Etrig_kev*anodecalib)-(aStep.Emesh_kev))/((aStep.Etrig_kev*anodecalib)+(aStep.Emesh_kev));
      delta_meshlength=((lengthenergy)-(aStep.Emesh_kev*meshcalib))/((lengthenergy)+(aStep.Emesh_kev*meshcalib));
      delta_anodelength=((aStep.Etrig_kev*anodecalib)-(lengthenergy))/((aStep.Etrig_kev*anodecalib)+(lengthenergy));
      delta_ccdlength=((aStep.Etrack_kev*weight)-(lengthenergy))/((aStep.Etrack_kev*weight)+(lengthenergy));

      if(aStep.EventNum-aStep.LastSpark<7)
	continue;

      if(aStep.Track_range_pix*MMPERPIXEL<=0.)
	continue;

      if(aStep.Track_maxpixel_ccdadu<72 ||
				     aStep.Track_maxpixel_ccdadu>200)
	continue;

      if(aStep.Etrig_kev*anodecalib<400)
	continue;

      if(aStep.Mesh_max_v>0.38)
	continue;

      if(aStep.Mesh_peak_v>0.38)
	continue;

      hist_ccd_vs_length->Fill(aStep.Etrack_kev*weight, lengthenergy);
      
    }

  Double_t fitf(Double_t *x, Double_t *par)
  {
    Double_t fitval = par[0]*x[0];
    return fitval;
  }   

  TF1 *func1 = new TF1("pol2", "pol2", 200, 10000);

  new TCanvas;
  hist_ccd_vs_length->SetTitle("CCD Energy Calibration");
  hist_ccd_vs_length->SetXTitle("Raw CCD Energy/keV");
  hist_ccd_vs_length->SetYTitle("Length Energy/keV");
  hist_ccd_vs_length->SetMarkerStyle(6);
  hist_ccd_vs_length->Draw("COLZ");
  hist_ccd_vs_length->Fit(func1, "r n");
  func1->SetLineColor(kRed);
  func1->Draw("SAME");

  cout<<"CCD Function Calibration : "<<(func1->GetParameter(2))<<"rawccd^2 + "<<(func1->GetParameter(1))<<"rawccd + "<<(func1->GetParameter(0))<<endl;

}

double fitPoly12Percent(double E)
{
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;                                                    
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;          

}    
