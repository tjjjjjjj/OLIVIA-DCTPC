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

//This macro is used to calibrate energy on a sequence-by-sequence basis

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
  
TH1D *hist_energy1=new TH1D("", "", 50, 0, 25000);
TH1D *hist_energy3=new TH1D("", "", 50, 0, 25000);

TH1D *hist_test=new TH1D("", "", 100, .0, 1.3);


double Etrig_kev_calib,Emesh_kev_calib;

  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;

int seq=0;
int runnum=0;
int color=1;
for (int event = 0; event<nentries; event++)
{

aStep.GetEntry(event); 

//sequence-by-sequence energy calibration
Etrig_kev_calib = aStep.Etrig_kev * ENERGY_CALIB_wf[aStep.SequenceNum]; 
Emesh_kev_calib = aStep.Emesh_kev * ENERGY_CALIB_wf[aStep.SequenceNum];


if(((aStep.SequenceNum!=seq) && event!=0)||event==nentries-1)
{
hist_energy2=(TH1D*)hist_energy1->Clone();
hist_energy2->SetLineColor(color);
hist_energy2->SetTitle("Event energies for each sequence (waveform only)");
hist_energy2->SetXTitle("Mesh energy (keV_{ee})");
hist_energy2->SetYTitle("counts/bin");
hist_energy2->Scale(1./EXPOSURE_SEQ[seq]);
hist_energy2->SetMaximum(0.1);
hist_energy2->Draw("SAME");
// gPad->GetCanvas()->SetLogy();
hist_energy1->Delete();
TH1D *hist_energy1=new TH1D("E", "E", 50,0, 25000);
}

hist_test->Fill(aStep.Mesh_base_v);


// if(aStep.SequenceNum!=seq)
// new TCanvas;

// if(aStep.SequenceNum==12)
// continue;

// if(aStep.Ntrack!=1)
// continue;
// 
// if(aStep.Anode_base_v>1.03)
// continue;
// 
// if(aStep.Mesh_base_v>0.2)
// continue;
// 
// if(aStep.Mesh_peak_v/Emesh_kev_calib>0.00015)
// continue;

// color=aStep.SequenceNum-6;
color=aStep.SequenceNum-11;
runnum=aStep.RunNum;

// double SRIMpar[3] = {0.14,0.9110,1.3};
// double Vdrift=1.15e06;
// double convertSample2sec = 4.e-09;//second per sample (250 MS/s)
// double Lvertical = aStep.Anode_R0time_samp*convertSample2sec*Vdrift;
// double Ltrack = sqrt((aStep.Track_range_pix*MMPERPIXEL/10.)*(aStep.Track_range_pix*MMPERPIXEL/10.)+(Lvertical)*(Lvertical));
// double deltaPoly = SRIMpar[1]*SRIMpar[1]-4*(SRIMpar[0]*(SRIMpar[2]-Ltrack));
// double Elength = (-1*SRIMpar[1]+sqrt(deltaPoly))/(2*SRIMpar[0]);



if(aStep.SequenceNum>11)
hist_energy1->Fill(Emesh_kev_calib);	
		
cout<<aStep.RunNum<<" "<<Emesh_kev_calib<<" "<<seq<<endl;		
		
		
		
seq=aStep.SequenceNum;	   
}

new TCanvas;
seq=0;


for (int event = 0; event<nentries; event++)
{

aStep.GetEntry(event);

//sequence-by-sequence energy calibration
Etrig_kev_calib = aStep.Etrig_kev * ENERGY_CALIB_wf[aStep.SequenceNum]; 
Emesh_kev_calib = aStep.Emesh_kev * ENERGY_CALIB_wf[aStep.SequenceNum];


if(((aStep.SequenceNum!=seq) && event!=0)||event==nentries-1)
{
hist_energy4=(TH1D*)hist_energy3->Clone();
hist_energy4->SetLineColor(color);
hist_energy4->SetTitle("Event energies for each sequence (waveform only)");
hist_energy4->SetXTitle("Anode energy (keV_{ee})");
hist_energy4->SetYTitle("counts/bin");
hist_energy4->Scale(1./EXPOSURE_SEQ[seq]);
hist_energy4->SetMaximum(0.1);
hist_energy4->Draw("SAME");
// gPad->GetCanvas()->SetLogy();
hist_energy3->Delete();
TH1D *hist_energy3=new TH1D("E", "E", 50,0, 25000);
}

// if(aStep.SequenceNum!=seq)
// new TCanvas;

// if(aStep.SequenceNum==12)
// continue;



// if(aStep.Ntrack!=1)
// continue;

if(aStep.Anode_base_v>1.03)
continue;

if(aStep.Mesh_base_v>0.2)
continue;

if(aStep.Mesh_peak_v/Emesh_kev_calib>0.00015)
continue;


//color=aStep.SequenceNum-6;
color=aStep.SequenceNum-11;
runnum=aStep.RunNum;
if(aStep.SequenceNum>11)
hist_energy3->Fill(Etrig_kev_calib);		
seq=aStep.SequenceNum;	   
}



// new TCanvas;
// hist_test->Draw();

}
