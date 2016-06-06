#define DCTPC_runtree_cxx
#include "../input_files/DCTPC_runtree.h"
#include "../input_files/Calibration_constants_LittleDCTPC.h";
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TRandom3.h"
#include "TGraph.h"
#include "TMath.h"
#include "TArc.h"
#include "TVirtualFitter.h"
#include <iostream>
#include <fstream>
TGraph *gr;
using namespace std;

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

  ofstream run_exposure;
  run_exposure.open("Run_Exposure.txt");
  
  TH2D *hist_time=new TH2D("Exposure Time","Exposure Time",10000, 1374000000 ,1386000000,1000,0,1000);//this gives a plot for the amount of exposure time i.e. the amount of time we are actually recording the events for against the time it was measured for
  
  TFile *outtree = new TFile("$LittleDCTPC_physics_input");
  TTree *dctreepc = (TTree*)outtree->Get("dctpc_runinfo");
  
  DCTPC_runtree aStep(dctreepc);
  
  Long64_t nentries = dctreepc->GetEntries();
  cout << "Number of Entries: " << nentries << endl;
  
  int exposure2=0;
  int exposure3=0;
  int exposure4=0;
  int exposure5=0;
  int exposure6=0;
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
  int exposure17=0;
  int exposure18=0;

  for (int event = 0; event<nentries; event++)
    {
      aStep.GetEntry(event);

      run_exposure<<aStep.Exposure_sec<<" , ";
      
      if(aStep.RunNum>=9535 && aStep.RunNum<=9906)
       	exposure2+=aStep.Exposure_sec;

      if(aStep.RunNum>=9907 && aStep.RunNum<=10001)
	exposure3+=aStep.Exposure_sec;

      if(aStep.RunNum>=10004 && aStep.RunNum<=10273)
	exposure4+=aStep.Exposure_sec;
      
      if(aStep.RunNum>=10274 && aStep.RunNum<=10671)
	exposure5+=aStep.Exposure_sec;

      if(aStep.RunNum>=10672 && aStep.RunNum<=11064)
	exposure6+=aStep.Exposure_sec;

      if(aStep.RunNum>=11065 && aStep.RunNum<=11161)
	exposure7+=aStep.Exposure_sec;

      if(aStep.RunNum>=11162 && aStep.RunNum<=11435)
	exposure8+=aStep.Exposure_sec;

      if(aStep.RunNum>=11436 && aStep.RunNum<=11788)
	exposure9+=aStep.Exposure_sec;

      if(aStep.RunNum>=11789 && aStep.RunNum<=12175)
	exposure10+=aStep.Exposure_sec;

      if(aStep.RunNum>=12176 && aStep.RunNum<=12561)
	exposure11+=aStep.Exposure_sec;

      if(aStep.RunNum>=12562 && aStep.RunNum<=12947)
	exposure12+=aStep.Exposure_sec;

      if(aStep.RunNum>=12948 && aStep.RunNum<=13340)
	exposure13+=aStep.Exposure_sec;

      if(aStep.RunNum>=13341 && aStep.RunNum<=13723)
	exposure14+=aStep.Exposure_sec;

      if(aStep.RunNum>=13724 && aStep.RunNum<=14120)
	exposure15+=aStep.Exposure_sec;

      if(aStep.RunNum>=14121 && aStep.RunNum<=14256)
	exposure16+=aStep.Exposure_sec;

      if(aStep.RunNum>=14262 && aStep.RunNum<=14374)
	exposure17+=aStep.Exposure_sec;

      if(aStep.RunNum>=14375 && aStep.RunNum<=15579)
	exposure18+=aStep.Exposure_sec;

      hist_time->Fill(aStep.Time_startofrun_sec,aStep.Exposure_sec); 
      
    }
  
  cout<<"Sequence 2 has an exposure of "<<exposure2<<" seconds"<<endl; 
  cout<<"Sequence 3 has an exposure of "<<exposure3<<" seconds"<<endl;
  cout<<"Sequence 4 has an exposure of "<<exposure4<<" seconds"<<endl;
  cout<<"Sequence 5 has an exposure of "<<exposure5<<" seconds"<<endl;
  cout<<"Sequence 6 has an exposure of "<<exposure6<<" seconds"<<endl;
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
  cout<<"Sequence 17 has an exposure of "<<exposure17<<" seconds"<<endl;
  cout<<"Sequence 18 has an exposure of "<<exposure18<<" seconds"<<endl;
  
  new TCanvas;
  hist_time->SetMarkerStyle(7);
  hist_time->SetTitle("Exposure Time");
  hist_time->SetXTitle("Unix time");
  hist_time->SetYTitle("Exposure / sec");
  hist_time->Draw();

  run_exposure.close();
  
}
