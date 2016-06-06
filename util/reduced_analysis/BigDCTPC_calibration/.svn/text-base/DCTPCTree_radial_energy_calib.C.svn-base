#define DCTPCTree_cxx
#include "../input_files/DCTPCTree.h"
#include "../input_files/Calibration_constants_BigDCTPC.h";
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


//This macro is used to study (and calibrate) the radial dependence of the light signal. Tracks that originate at the center of the TPC are brighter than those on the side and this needs to be accounted for in energy reconstruction.
void DCTPCTree::Loop()
{

TStopwatch timer;
  gROOT->SetStyle("Plain");
  gStyle->SetEndErrorSize(3);
  gStyle->SetPalette(1,0);
  gStyle->SetLineWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(kTRUE);
  TH1::AddDirectory(false);
  
  TFile *outtree = new TFile("$BigDCTPC_calibration_input");
  TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
  DCTPCTree aStep(dctreepc);
  
  TH2D *hist_radialE=new TH2D("","",30,0,600,30,0,2.);
  TH2D *hist_energy=new TH2D("E", "E", 100, 0, 1000, 100, 0, 1000);
  TH2D *hist_energy2=new TH2D("E", "E", 100, 0, 1000, 100, 0, 1000);
  
double recalibmesh = 1.;
double recalibCCD = 1.;
double recalibanode = 1.;
double rstart=0.;
double rend=0.;

  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;
 
  for (int event = 0; event<nentries; event++)
    {
      aStep.GetEntry(event); 
      
      if(event%10000==0)
	  cout<<((double)event/(double)nentries)*100.<<"%"<<endl;	

double rr=(pow(aStep.Track_x_pix,2)+pow(aStep.Track_y_pix,2))*pow(MMPERPIXEL/10.,2);

 //these cuts are meant to isolate short, fully contained tracks not originating from the rings. We want short tracks so that the radial position is well defined.
     if (
     aStep.Edge==0  && 
     aStep.Ntrig <=2 &&
TMath::Abs((aStep.Etrig_kev-aStep.Etrack_kev)/(aStep.Etrig_kev+aStep.Etrack_kev))<0.4 &&
TMath::Abs((aStep.Etrig_kev-aStep.Emesh_kev)/(aStep.Etrig_kev+aStep.Emesh_kev))<0.1 &&
     aStep.Track_range_pix<80.
     )	
	{
 

	 hist_radialE->Fill(rr,(aStep.Etrack_kev)/(aStep.Etrig_kev));
hist_energy->Fill(aStep.Etrack_kev,aStep.Etrig_kev);  
hist_energy2->Fill(aStep.Etrig_kev,aStep.Emesh_kev);
	 
	 }
    
     }
 
new TCanvas;
hist_radialE->SetXTitle("radius (pixels)");
hist_radialE->SetYTitle("E_{CCD}/E_{anode}");
hist_radialE->Draw("COLZ");  
  
TProfile *prof = hist_radialE->ProfileX();
new TCanvas;

prof->SetYTitle("E_{CCD}/E_{anode}");
prof->Draw();  


//perform a simple fit that can be used for calibrating energy as a function fo radius  
TF1 *f1 = new TF1("f1","[0]/([1] +x)^[2]",0,440);
prof->Fit("f1","R");
  
cout<<"Param 0: "<<f1->GetParameter(0)<<endl;
cout<<"Param 1: "<<f1->GetParameter(1)<<endl;
cout<<"Param 2: "<<f1->GetParameter(2)<<endl;

//6.23665/pow(141.035+r,0.360402);

TF1 *f2 = new TF1("f2","7.56874/(146.405+x)^0.388667",0,440); 
  
f2->Draw("SAME");
new TCanvas;
hist_energy->Draw("COLZ"); 

new TCanvas;
hist_energy2->Draw("COLZ");   
  
}
