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

//This macro is used for determining the efficiency of the dd source (and Po210) calibration runs

void DCTPCTree::Loop()
{

TStopwatch timer;
gROOT->SetStyle("Plain");
gStyle->SetEndErrorSize(3);
gStyle->SetPalette(1,0);
gStyle->SetLineWidth(2);
gStyle->SetHistLineWidth(2);
gStyle->SetOptStat(kFALSE);
// gStyle->SetOptFit(kFALSE);
TH1::AddDirectory(false);
  
TFile *outtree = new TFile("$BigDCTPC_calibration_input");
TTree *dctreepc = (TTree*)outtree->Get("dctpc_eventinfo");
DCTPCTree aStep(dctreepc);

TH2D *hist_energy=new TH2D("", "", 50, 0, 10000, 50, 0, 10000);

TH1D *hist_energyeff=new TH1D("", "",40,0,8000);
TH1D *hist_energyeff_wf=new TH1D("", "",40,0,8000);
TH1D *hist_mesh=new TH1D("", "",40,0,8000);

double rstart=0.;
double rend=0.;
int points=0;
int points2=0;
double x[1000000],y[1000000],x2[1000000],y2[1000000];
double Etrack_kev_calib,Etrig_kev_calib,Emesh_kev_calib,TrackX_pix_calib,TrackY_pix_calib,radius,TrackX_pix_calib_start,TrackY_pix_calib_start,TrackX_pix_calib_end,TrackY_pix_calib_end;
//create event loop
  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;
 
for (int event = 0; event<nentries; event++)
{
      aStep.GetEntry(event); 
      
      if(event%100000==0)
	  cout<<((double)event/(double)nentries)*100.<<"%"<<endl;	

//calibration of variables 

//set-by-set position calibration
TrackX_pix_calib = aStep.Track_x_pix - POSITION_X_CALIB[aStep.SetNum];
TrackY_pix_calib = aStep.Track_y_pix - POSITION_Y_CALIB[aStep.SetNum];
 
TrackX_pix_calib_start = aStep.Track_x_start_pix - POSITION_X_CALIB[aStep.SetNum];
TrackY_pix_calib_start = aStep.Track_y_start_pix - POSITION_Y_CALIB[aStep.SetNum];
TrackX_pix_calib_end = aStep.Track_x_end_pix - POSITION_X_CALIB[aStep.SetNum];
TrackY_pix_calib_end = aStep.Track_y_end_pix - POSITION_Y_CALIB[aStep.SetNum];

//sequence-by-sequence energy calibration
Etrack_kev_calib = aStep.Etrack_kev * ENERGY_CALIB_ccd[aStep.SequenceNum];
Etrig_kev_calib = aStep.Etrig_kev * ENERGY_CALIB_wf[aStep.SequenceNum]; 
Emesh_kev_calib = aStep.Emesh_kev * ENERGY_CALIB_wf[aStep.SequenceNum];

//radial-based energy calibration
radius=(pow(TrackX_pix_calib,2)+pow(TrackY_pix_calib,2))*pow(MMPERPIXEL/10.,2);  
Etrack_kev_calib*=(pow(radius+RADIAL_PARAM1_CALIB,RADIAL_PARAM2_CALIB)/RADIAL_PARAM0_CALIB);

rstart=sqrt(pow(TrackX_pix_calib_start,2)+pow(TrackY_pix_calib_start,2));
rend=sqrt(pow(TrackX_pix_calib_end,2)+pow(TrackY_pix_calib_end,2));


// hist_energy->Fill(Emesh_kev_calib,Etrack_kev_calib);

//the efficiency reported does not (yet) include the fact that the CCD cannot see the whole field cage.

	 if(Etrig_kev_calib>0.&&Emesh_kev_calib>0.&&
	 Etrig_kev_calib<10000.&&Emesh_kev_calib<10000. 	 
&& TMath::Abs((Etrig_kev_calib-Emesh_kev_calib)/(Etrig_kev_calib+Emesh_kev_calib))<0.2
&& aStep.Ntrig<=2
	 )
	{

	 hist_energyeff_wf->Fill(Emesh_kev_calib);

	}


	 if(Etrig_kev_calib>0.&&Emesh_kev_calib>0.&&Etrack_kev_calib>0.&&
	 Etrig_kev_calib<10000.&&Emesh_kev_calib<10000.&&Etrack_kev_calib<10000. 
	 
&& TMath::Abs((Etrig_kev_calib-Emesh_kev_calib)/(Etrig_kev_calib+Emesh_kev_calib))<0.2
&& TMath::Abs((Etrack_kev_calib-Emesh_kev_calib)/(Etrack_kev_calib+Emesh_kev_calib))<0.4
&& Etrack_kev_calib>0.
&& aStep.Ntrig<=2
	 )
	{
	
   hist_mesh->Fill(Emesh_kev_calib);
   
	}
   
}


for(int i=1;i<=hist_energyeff->GetNbinsX();i++)
{
hist_energyeff->SetBinContent(i,hist_mesh->GetBinContent(i)/hist_energyeff_wf->GetBinContent(i));
hist_energyeff->SetBinError(i,hist_energyeff->GetBinContent(i)*sqrt((1./hist_mesh->GetBinContent(i))+(1./hist_energyeff_wf->GetBinContent(i))));
}


//account for the fact that the anode (20 cm radius) is bigger than the imaging area 
hist_energyeff->Scale((3.14159*pow(20.,2))/POSITION_IMAGING_AREA[0]);


new TCanvas;
hist_energyeff->SetMinimum(0.);

hist_energyeff->Draw();
hist_energyeff->SetTitle("Detection efficiency");
hist_energyeff->SetXTitle("Mesh energy (keVee)");
hist_energyeff->SetYTitle("Efficiency fraction");

new TCanvas;
hist_energyeff_wf->Draw();
hist_mesh->SetLineColor(2);
hist_mesh->Draw("SAME");



}
