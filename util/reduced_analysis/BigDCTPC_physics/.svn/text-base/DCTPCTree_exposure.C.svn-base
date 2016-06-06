#define DCTPC_runtree_cxx
#include "../input_files/DCTPC_runtree.h"
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

//This macro is used to find the exposure in seconds for every sequence

void DCTPC_runtree::Loop()
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

TH2D *hist_time=new TH2D("","",10000,1391000000,1410000000,1000,0,1000);
TH2D *hist_time2=new TH2D("","",10000,900,15200,1000,0,1000);
TFile *outtree = new TFile("$BigDCTPC_physics_input");
TTree *dctreepc = (TTree*)outtree->Get("dctpc_runinfo");

DCTPC_runtree aStep(dctreepc);

Long64_t nentries = dctreepc->GetEntries();
cout << "Number of Entries: " << nentries << endl;

int exposure7=0;
int exposure8=0;
int exposure9=0;
int exposure10=0;
int exposure11=0;
int exposure12=0;
int exposure13=0;
int exposure14=0;
int exposure15=0;
int exposure16=0;

for (int event = 0; event<nentries; event++)
{
aStep.GetEntry(event);

 
 if(aStep.RunNum>=955 && aStep.RunNum<=1645)
 exposure7+=aStep.Exposure_sec;
 
 if(aStep.RunNum>=1646 && aStep.RunNum<=2205)
 exposure8+=aStep.Exposure_sec; 
 
 if(aStep.RunNum>=2206 && aStep.RunNum<=2648)
 exposure9+=aStep.Exposure_sec;

 if(aStep.RunNum>=2650 && aStep.RunNum<=3868)
 exposure10+=aStep.Exposure_sec;
 
 if(aStep.RunNum>=3876 & aStep.RunNum<=5789)
 exposure11+=aStep.Exposure_sec;   
 
 if(aStep.RunNum>=5792&&aStep.RunNum<=7520)
 exposure12+=aStep.Exposure_sec;
 
 if(aStep.RunNum>=7530&&aStep.RunNum<=9143)
 exposure13+=aStep.Exposure_sec;
 
 if(aStep.RunNum>=9147&&aStep.RunNum<=12833)
 exposure14+=aStep.Exposure_sec;
 
 if(aStep.RunNum>=12837&&aStep.RunNum<=13311)
 exposure15+=aStep.Exposure_sec;
 
 if(aStep.RunNum>=13315)
 exposure16+=aStep.Exposure_sec;
 
hist_time->Fill(aStep.Time_startofrun_sec,aStep.Exposure_sec); 
hist_time2->Fill(aStep.RunNum,aStep.Exposure_sec); 
}

cout<<"Sequence 7 has an exposure of "<<exposure7<<" seconds"<<endl; 
cout<<"Sequence 8 has an exposure of "<<exposure8<<" seconds"<<endl;
cout<<"Sequence 9 has an exposure of "<<exposure9<<" seconds"<<endl;
cout<<"Sequence 10 has an exposure of "<<exposure10<<" seconds"<<endl;
cout<<"Sequence 11 has an exposure of "<<exposure11<<" seconds"<<endl;
cout<<"Sequence 12 has an exposure of "<<exposure12<<" seconds"<<endl;
cout<<"Sequence 13 has an exposure of "<<exposure13<<" seconds"<<endl;
cout<<"Sequence 14 has an exposure of "<<exposure14<<" seconds"<<endl;
cout<<"Sequence 15 has an exposure of "<<exposure15<<" seconds"<<endl;
cout<<"Sequence 16 has an exposure of "<<exposure16<<" seconds"<<endl;

new TCanvas;
hist_time->SetMarkerStyle(7);
hist_time->SetTitle("Exposure for each run");
hist_time->SetXTitle("Unix time");
hist_time->SetYTitle("Exposure (seconds)");
hist_time->Draw();


new TCanvas;
hist_time2->SetMarkerStyle(7);
hist_time2->SetTitle("Exposure for each run");
hist_time2->SetXTitle("Unix time");
hist_time2->SetYTitle("Exposure (seconds)");
hist_time2->Draw();

}
