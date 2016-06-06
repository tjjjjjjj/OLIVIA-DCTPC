#define DCTPCTree_cxx
#include "../input_files/DCTPCTree.h"
#include "../input_files/Calibration_constants_LittleDCTPC.h";
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//This macro is used to study (and calibrate) the radial dependence of the light signal. Tracks that originate at the center of the TPC are brighter than those on the side and this needs to be accounted for in energy reconstruction.
//The neutron calibration is very useful for this purpose as the tracks are quite short, which makes the radial position of the track well defined.
void DCTPCTree::Loop()
{

TStopwatch timer;
  gROOT->SetStyle("Plain");
  gStyle->SetEndErrorSize(3);
  gStyle->SetPalette(1,0);
  gStyle->SetLineWidth(2);
  gStyle->SetHistLineWidth(2);
//   gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(kTRUE);
  TH1::AddDirectory(false);
  
  TFile *outtree = new TFile("$LittleDCTPC_neutron_calibration_input");
  TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
  DCTPCTree aStep(dctreepc);
  
  TH2D *hist_radialE=new TH2D("","",30,0,150,30,0,2.0);
  TH2D *hist_energy=new TH2D("E", "E", 100, 0, 1000, 100, 0, 1000);

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
     if (aStep.BurnIn==0 && 
     aStep.Edge==0  && 
     (aStep.Etrack_kev-aStep.Etrig_kev)<50. && 
     aStep.Track_range_pix<80. )	
	{
	
     //this cut removes the events that are associated with a very small waveform, probably due to noise
	 if(abs((aStep.Etrig_kev)-(aStep.Etrack_kev))>1000. && (aStep.Etrig_kev)<2000.)
	 continue;

	 hist_radialE->Fill(rr,(aStep.Etrack_kev)/(aStep.Etrig_kev));	 	
	 hist_energy->Fill(aStep.Etrack_kev,aStep.Etrig_kev);	
	 	  
	 }
    
     }
 
new TCanvas;
hist_radialE->Draw("COLZ");  
  
TProfile *prof = hist_radialE->ProfileX();
new TCanvas;
prof->Draw();  


//perform a simple fit that can be used for calibrating energy as a function of radius  
TF1 *f1 = new TF1("f1","[0]/([1] +x)^[2]",0,500);
prof->Fit("f1");
  
cout<<"Param 0: "<<f1->GetParameter(0)<<endl;
cout<<"Param 1: "<<f1->GetParameter(1)<<endl;
cout<<"Param 2: "<<f1->GetParameter(2)<<endl;

//this is the original fit used in the DCTPC NIM paper. It is interesting to compare the two fits. 
TF1 *f2 = new TF1("f2","1./(pow(1+ x/(16.123*16.123),3.5857/2 ))",0,500);
f2->SetLineColor(2);
f2->Draw("SAME");

new TCanvas; 
f2->Draw();
f1->Draw("SAME");

new TCanvas;
hist_energy->Draw("COLZ");
  
}
