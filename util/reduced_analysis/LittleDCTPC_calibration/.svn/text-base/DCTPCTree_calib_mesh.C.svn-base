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
  
  TFile *outtree = new TFile("../input_files/allie.root");
  TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
  DCTPCTree aStep(dctreepc);


   TH2D *hist_emesh=new TH2D("Emesh", "Emesh", 1000, 0, 10000, 1000, 0, 10000);
   TH1D *hist_delta=new TH1D("delta","delta",50,-1.,1.);


  
//create event loop
  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;

  for (int event = 0; event<nentries; event++)
    {
      aStep.GetEntry(event);
      if(event%10000==0)
      cout<<((double)event/(double)nentries)*100.<<"%"<<endl;

      double r2 = (aStep.Track_x_pix*aStep.Track_x_pix+aStep.Track_y_pix*aStep.Track_y_pix)*pow(MMPERPIXEL/10.,2);
       double weight=1./pow(RADIAL_PARAM0_CALIB/(RADIAL_PARAM1_CALIB+r2),RADIAL_PARAM2_CALIB);//radial calibration 
       double anodecalib=0.958675;

      double delta_anodeccd=((aStep.Etrig_kev*anodecalib)-(aStep.Etrack_kev*weight))/((aStep.Etrig_kev*anodecalib)+(aStep.Etrack_kev*weight));
      double delta_meshccd=((aStep.Emesh_kev)-(aStep.Etrack_kev*weight))/((aStep.Emesh_kev)+(aStep.Etrack_kev*weight));
      double delta_anodemesh=((aStep.Etrig_kev*anodecalib)-(aStep.Emesh_kev))/((aStep.Etrig_kev*anodecalib)+(aStep.Emesh_kev));      
      
            
      //if(aStep.EventNum-aStep.LastSpark<7)
      //continue;

      //if(aStep.Track_range_pix*MMPERPIXEL<=0.)
      //continue;

      //if(aStep.Track_maxpixel_ccdadu>200)
      //continue;

      //if(fabs(aStep.Mesh_R0time_samp)>0.5 ||
      // fabs(aStep.Mesh_R10time_samp)>0.5 ||
      // fabs(aStep.Mesh_R25time_samp)>0.5 ||
      // fabs(aStep.Mesh_R50time_samp)>0.5 ||
      // fabs(aStep.Mesh_R75time_samp)>0.5 ||
      // fabs(aStep.Mesh_R90time_samp)>0.5)
      //continue;

      // if(fabs(aStep.Mesh_F0time_samp)>0.5 ||
      // fabs(aStep.Mesh_F10time_samp)>0.5 ||
      // fabs(aStep.Mesh_F25time_samp)>0.5 ||
      // fabs(aStep.Mesh_F50time_samp)>0.5 ||
      // fabs(aStep.Mesh_F75time_samp)>0.5 ||
      // fabs(aStep.Mesh_F90time_samp)>0.5)
      //continue;

      //if(aStep.Etrig_kev<1250 &&
      //TMath::Abs(delta_anodeccd)>0.4)
      //continue;

      // if(TMath::Abs(delta_meshccd)>0.4)
      //continue;

      if(aStep.Etrig_kev<1250 &&
			 &&aStep.Mesh_max_v<0.39 &&
					    TMath::Abs(delta_anodemesh)>0.4)
      continue;


	  hist_delta->Fill((aStep.Etrig_kev*(anodecalib)-(aStep.Emesh_kev))/(aStep.Etrig_kev*(anodecalib)+(aStep.Emesh_kev)));

	  hist_emesh->Fill(aStep.Emesh_kev, aStep.Etrig_kev*(anodecalib));

    }

  Double_t fitf(Double_t *x, Double_t *par)
  {
    Double_t fitval = par[0]*x[0];
    return fitval;
  }   

  // TF1 *func = new TF1("fitf",fitf, 400, 1900,1);
  //func->SetParameter(0,1);
  //hist_emesh->Fit(func, "r n");
  TF1 *func = new TF1("pol1", "pol1", 700, 1800);

  new TCanvas;
  hist_emesh->SetTitle("Mesh Energy Calibration");
  hist_emesh->SetXTitle("Raw Mesh Energy/keV");
  hist_emesh->SetYTitle("Calibrated Anode Energy/keV");
  hist_emesh->SetMarkerStyle(6);
  hist_emesh->Draw("COLZ");
  hist_emesh->Fit(func, "r n");
  func->SetLineColor(kRed);
  func->Draw("SAME");

  new TCanvas;
  hist_delta->SetTitle("Anode Energy Minus Mesh Energy");
  hist_delta->SetXTitle("Calibrated Anode Energy/keV Minus Raw Mesh Energy/keV");
  hist_delta->SetYTitle("Counts");
   hist_delta->Draw();
  hist_delta->Fit("gaus");
  gaus->SetLineColor(kBlue);
  gaus->Draw("SAME");


  cout<<"Anode Energy - Mesh Energy= "<<gaus->GetParameter(1)<<endl;
  cout<<"Calibrated Anode Energy / Raw Mesh Energy: "<<func->GetParameter(1)<<endl;
  cout<<"Mesh Energy Calibration Ratio: "<<func->GetParameter(1)<<endl;
  cout<<"Anode Mesh Linear Fit Parameter 0 : "<<func->GetParameter(0)<<endl;
  cout<<"Anode Mesh Linear Fit Parameter 1 : "<<func->GetParameter(1)<<endl;
  
}

double fitPoly12Percent(double E)
{
  return sqrt(0.9741*0.9741*E*E+113.835*113.835)-62.675;
  //return -11.50+39.02*E+0.03459*E*E;                                                    
  //return 469 + 10*E + 0.548 *E*E -3.39e-3*E*E*E+7.5e-6*E*E*E*E;          

}    
