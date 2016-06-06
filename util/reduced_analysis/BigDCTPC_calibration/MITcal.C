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
TGraph *gr;


void DCTPCTree::Loop()
{

TStopwatch timer;
gROOT->SetStyle("Plain");
gStyle->SetEndErrorSize(3);
gStyle->SetPalette(1,0);
gStyle->SetLineWidth(2);
gStyle->SetHistLineWidth(2);
gStyle->SetOptStat(1);
gStyle->SetOptFit(1);
TH1::AddDirectory(false);
  
TFile *outtree = new TFile("$BigDCTPC_calibration_input");
TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
DCTPCTree aStep(dctreepc);

TH1D *hist_energy_ccd=new TH1D("CCD energy","CCD energy",100,10,15000);
TH1D *hist_energy_trig=new TH1D("Anode energy","Anode energy",100,10,15000);
TH1D *hist_energy_mesh=new TH1D("Mesh energy","Mesh energy",100,10,15000);

TH2D *hist_pos=new TH2D("pos","pos",256,-612,612,256,-612,612);

double recalibCCD[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0}; 	//{0.7663,0.8900,0.8887,0.8644,0.7078,0.968,0.9719};
double recalibanode[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};	//{0.7710,0.8861,0.8849,0.8570,0.7150,0.8401,0.8291};
double recalibmesh[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};	//{0.769,0.8853,0.8810,0.856,0.7134,0.8392,0.8284};
int j;

double rstart=0.;
double rend=0.;
double ralpha=0.;


//create event loop
  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;
 
for (int event = 0; event<nentries; event++)
{
    aStep.GetEntry(event); 
      
    if(event%100000==0)
    {
		cout<<((double)event/(double)nentries)*100.<<"%"<<endl;	
	}
	
	if      (aStep.RunNum>=889 && aStep.RunNum<=891){j=0;}
	else if (aStep.RunNum>=893 && aStep.RunNum<=896){j=1;}
	else if (aStep.RunNum>=897 && aStep.RunNum<=904){j=2;}
	else if (aStep.RunNum>=910 && aStep.RunNum<=917){j=3;}
	else if (aStep.RunNum>=919 && aStep.RunNum<=925){j=4;}	
	else if (aStep.RunNum>=926 && aStep.RunNum<=938){j=5;}	
	else if (aStep.RunNum>=939 && aStep.RunNum<=954){j=6;}	
	else {continue;}

	rstart=sqrt(pow(aStep.Track_x_start_pix,2)+pow(aStep.Track_y_start_pix,2));	
	ralpha=sqrt(pow(aStep.Track_x_start_pix+313.125,2)+pow(aStep.Track_y_start_pix+55.625,2));
 
      //for calculating everything if it passes.
     if ( aStep.Edge==0   && (aStep.Etrack_kev-aStep.Etrig_kev)<500000. && ralpha<15) //15
	 {
	 	hist_pos->Fill(aStep.Track_x_start_pix,aStep.Track_y_start_pix);
	  
	 	hist_energy_ccd->Fill(aStep.Etrack_kev*recalibCCD[j]); 
	 	hist_energy_trig->Fill(aStep.Etrig_kev*recalibanode[j]); 
	 	hist_energy_mesh->Fill(aStep.Emesh_kev*recalibmesh[j]);
	 }
}

new TCanvas;
hist_pos->Draw("COLZ");

TF1 *f1 = new TF1("f1","gaus",3700,6900);

new TCanvas;
hist_energy_ccd->SetTitle("CCD Energy");
hist_energy_ccd->SetXTitle("keV");
hist_energy_ccd->SetYTitle("Counts");
hist_energy_ccd->Draw();
hist_energy_ccd->Fit("f1","R");

TF1 *f2 = new TF1("f2","gaus",4300,6300);

new TCanvas;
hist_energy_trig->SetTitle("Anode Energy");
hist_energy_trig->SetXTitle("keV");
hist_energy_trig->SetYTitle("Counts");
hist_energy_trig->Draw();
hist_energy_trig->Fit("f2","R");

TF1 *f3 = new TF1("f3","gaus",4300,6300);

new TCanvas;
hist_energy_mesh->SetTitle("Mesh Energy");
hist_energy_mesh->SetXTitle("keV");
hist_energy_mesh->SetYTitle("Counts");
hist_energy_mesh->Draw();
hist_energy_mesh->Fit("f3","R");

cout<<"*******************"<<endl;
cout<<" "<<endl;
cout<<"Sequence 1 CCD scale: "<<recalibCCD[0]<<" Anode scale: "<<recalibanode[0]<<" Mesh scale: "<<recalibmesh[0]<<endl;
cout<<"Sequence 2 CCD scale: "<<recalibCCD[1]<<" Anode scale: "<<recalibanode[1]<<" Mesh scale: "<<recalibmesh[1]<<endl;
cout<<"Sequence 3 CCD scale: "<<recalibCCD[2]<<" Anode scale: "<<recalibanode[2]<<" Mesh scale: "<<recalibmesh[2]<<endl;
cout<<"Sequence 4 CCD scale: "<<recalibCCD[3]<<" Anode scale: "<<recalibanode[3]<<" Mesh scale: "<<recalibmesh[3]<<endl;
cout<<"Sequence 5 CCD scale: "<<recalibCCD[4]<<" Anode scale: "<<recalibanode[4]<<" Mesh scale: "<<recalibmesh[4]<<endl;
cout<<"Sequence 6 CCD scale: "<<recalibCCD[5]<<" Anode scale: "<<recalibanode[5]<<" Mesh scale: "<<recalibmesh[5]<<endl;
cout<<"Sequence 7 CCD scale: "<<recalibCCD[6]<<" Anode scale: "<<recalibanode[6]<<" Mesh scale: "<<recalibmesh[6]<<endl;
cout<<" "<<endl;
cout<<"CCD fit sigma/mean: "<<f1->GetParameter(2)/f1->GetParameter(1)<<endl;
cout<<"Anode fit sigma/mean: "<<f2->GetParameter(2)/f2->GetParameter(1)<<endl;
cout<<"Mesh fit sigma/mean: "<<f3->GetParameter(2)/f3->GetParameter(1)<<endl;
}
