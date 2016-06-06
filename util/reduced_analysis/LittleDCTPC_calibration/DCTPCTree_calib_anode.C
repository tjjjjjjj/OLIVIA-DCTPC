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
  
  TFile *outtree = new TFile("$LittleDCTPC_alpha_calibration_input");
  TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
  DCTPCTree aStep(dctreepc);


  TH1D *hist_etrig=new TH1D("Eanode", "Eanode", 60, 0, 6000);
  TH1D *hist_anode_volt=new TH1D("Vanode", "Vanode", 100, 0, 1);


//create event loop
  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;

  for (int event = 0; event<nentries; event++)
    {
      aStep.GetEntry(event);
      if(event%10000==0)
      cout<<((double)event/(double)nentries)*100.<<"%"<<endl;

      if(event>0)

	{
	hist_etrig->Fill(aStep.Etrig_kev);
	hist_anode_volt->Fill(aStep.Anode_max_v);
	
	}

    }

  TF1 *f1 = new TF1("f1","gaus", 3900, 5100);
  hist_etrig->Fit("f1", "R N");
  TF1 *gaus1 = new TF1("gaus", "gaus", 0, 10);
  TF1 *gaus2 = new TF1("gaus", "gaus", 0, 10);


  cout<<"Mean Alpha Anode Energy Measured: "<<f1->GetParameter(1)<<endl;
  cout<<"Mean Alpha Anode Energy Source: 4400"<<endl;
  cout<<"Anode Energy Calibration Ratio: "<<4400/f1->GetParameter(1)<<endl;
  cout<<"Anode Calibration Gaus Fit Parameter 0 : "<<f1->GetParameter(0)<<endl;
  cout<<"Anode Calibration Gaus Fit Parameter 1 : "<<f1->GetParameter(1)<<endl;
  cout<<"Anode Calibration Gaus Fit Parameter 2 : "<<f1->GetParameter(2)<<endl;
//this is the constant you need to multiply all of the Etrig raw values by to get the calibrated Etrig values

  new TCanvas;
  hist_etrig->SetTitle("Anode Energy Calibration");
  hist_etrig->SetXTitle("Raw Anode Energy/keV");
  hist_etrig->SetYTitle("Counts");
  hist_etrig->Draw();
  f1->SetLineColor(kGreen);
  f1->Draw("SAME");

  new TCanvas;
  hist_anode_volt->SetTitle("Volt to keV Calibration");
  hist_anode_volt->SetXTitle("Voltage");
  hist_anode_volt->SetYTitle("Counts");
  hist_anode_volt->Draw();
  hist_anode_volt->Fit(gaus1, "r n");
  gaus1->SetLineColor(kRed);
  gaus1->Draw("SAME");
 
  cout<<"Mean Anode Maximum Voltage = "<<gaus1->GetParameter(1)<<endl;
  cout<<"Anode Voltage Calibration Ratio: "<<4400/gaus1->GetParameter(1)<<endl;
  
}
