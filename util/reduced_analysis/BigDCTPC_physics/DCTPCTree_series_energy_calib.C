#define DCTPCTree_cxx
#include "../input_files/DCTPCTree.h"
#include "../input_files/Calibration_constants_BigDCTPC.h";
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TRandom3.h"
#include "TGraph.h"
#include "TMath.h"
#include "TArc.h"
#include "TVirtualFitter.h"

//This macro is used to study the consistency between different series' (sets of 100 runs) within a sequence

void DCTPCTree::Loop()
{

TStopwatch timer;
gROOT->SetStyle("Plain");
gStyle->SetEndErrorSize(3);
gStyle->SetPalette(1,0);
gStyle->SetLineWidth(2);
gStyle->SetHistLineWidth(2);
gStyle->SetOptStat(kFALSE);
gStyle->SetOptFit(kFALSE);
TH1::AddDirectory(false);
  
TFile *outtree = new TFile("$BigDCTPC_physics_input");
TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
DCTPCTree aStep(dctreepc);
  
TH1D *hist_energy1=new TH1D("", "", 150, 0, 15000);

double Etrig_kev_calib,Emesh_kev_calib;

  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;

int ser=0;
int seq=0;
int runnum=0;
int color=1;
for (int event = 0; event<nentries; event++)
{

aStep.GetEntry(event); 

//sequence-by-sequence energy calibration
Etrig_kev_calib = aStep.Etrig_kev * ENERGY_CALIB_wf[aStep.SequenceNum]; 
Emesh_kev_calib = aStep.Emesh_kev * ENERGY_CALIB_wf[aStep.SequenceNum];

if(((aStep.SeriesNum!=ser) && event!=0)||event==nentries-1)
{
hist_energy2=(TH1D*)hist_energy1->Clone();
hist_energy2->SetLineColor(color);
hist_energy2->SetTitle("Event energies for each sequence (waveform only)");
hist_energy2->SetXTitle("Mesh energy (keV_{ee})");
hist_energy2->SetYTitle("counts/bin");
hist_energy2->Scale(1./EXPOSURE_SEQ[seq]);
hist_energy2->SetMaximum(0.1);
hist_energy2->Draw("SAME");
gPad->GetCanvas()->SetLogy();
hist_energy1->Delete();
TH1D *hist_energy1=new TH1D("E", "E", 150,0, 15000);
}

if(aStep.SequenceNum!=seq)
new TCanvas;

color=aStep.SequenceNum-6;
runnum=aStep.RunNum;
hist_energy1->Fill(Emesh_kev_calib);	
	
ser=aStep.SeriesNum;	
seq=aStep.SequenceNum;	   
}


}
