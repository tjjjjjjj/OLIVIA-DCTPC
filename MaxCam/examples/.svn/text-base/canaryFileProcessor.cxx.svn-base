// this is built off of looking at Jeremy's makeSimulations.cc

// headers
#include "TH2.h"
#include "TSystem.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCollection.h"

#include "../DmtpcDataset.hh"
#include "../MaxCamConfig.hh"
#include "../MaxCamImageTools.hh"
#include "../MaxCamChannel.hh"

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;



TH2F *bias=0;
void makeBiasFrame(TString fname) {
    if (bias) { delete bias; bias=0; }
    DmtpcDataset *d= new DmtpcDataset;
    d->openRootFile(fname);
	
    int nev=d->tree()->GetEntries();
    for (int i=0; i<nev; i++) {
        d->getEvent(i);
        if (!bias) bias=(TH2F*)d->event()->ccdData(0)->Clone("biasFrame1");
        else bias->Add( d->event()->ccdData(0) );
    }
    bias->Scale(1./nev);
    delete d; d=0;
}



// function prototypes
void fitsFileProcessor(int run) {

  
  // inputs
  int npix=1024;
  int nbins=512;

  cout << "inside fitsFileProcessor" << endl;

  cout << "Initializing..."<<endl;
  
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ////Set up data file                                                         //
  ///////////////////////////////////////////////////////////////////////////////
  DmtpcDataset *dataset = new DmtpcDataset();
  TString infile="data/canary_";
  if (run<10) infile+="000";
  else if (run<100) infile+="00";
  else if (run<1000) infile+="0";
  
  infile += run;
  
  TString outputfile=infile;
  outputfile += ".root";
  infile += ".fits";

  cout << "Input ...... " << infile << endl;
  cout << "Output ..... " << outputfile << endl;



  dataset->createRootFile( outputfile, "recreate");
  
  dataset->setComment("Calibration for canary; 0.5sec exposure; 720V, -1500V; emccd Luca");
  dataset->setLocation("44");
  dataset->setKeyword("KEYWORD");
  dataset->comment()->Write();
  dataset->keyword()->Write();
  dataset->location()->Write();

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ////Set up experimentConfig values. (copied from DmtpcRun.cc::L44-74)        //
  ///////////////////////////////////////////////////////////////////////////////
  cout << "Set up experimentConfig"<<endl;
  
  MaxCamChannel* driftHV   = new MaxCamChannel("driftHV",    "Drift Voltage", -1,  -1);
  driftHV->currentValue=driftHV->setValue=1.5;

  MaxCamChannel* anodeHV   = new MaxCamChannel("anodeHV",    "Anode Voltage", -1,  -1);
  anodeHV->currentValue = anodeHV->setValue = 0.68;

  MaxCamChannel* pressure = new MaxCamChannel("pressure",   "Gas pressure",  -1, -1);    
  pressure->currentValue = pressure->setValue = 75;

  cout << "Set up MaxCamConfig" << endl;
  
  MaxCamConfig *ccdconfig = new MaxCamConfig("ccdConfig","CCD Camera configuration");
  ccdconfig->cameraID = 1;
  ccdconfig->row_width = npix;
  ccdconfig->img_rows = npix;
  ccdconfig->hbin = npix / nbins;
  ccdconfig->vbin = npix / nbins;
  ccdconfig->ul_x = 0;
  ccdconfig->ul_y = 0;
  ccdconfig->lr_x = npix;
  ccdconfig->lr_y = npix;
  ccdconfig->CCDTemp = -20;
  ccdconfig->CCDTempSet = -20;
  ccdconfig->exposureTime = 1000;
  ccdconfig->frameType = 0;
  ccdconfig->nFlushes = -1;
  ccdconfig->bitDepth = 65535;
  ccdconfig->daqTime = -1;
  
  if (!bias)	{
	  bias = new TH2F("biasFrame1", "biasFrame1",  nbins,0,nbins,   nbins,0,nbins);
	  for (int i = 1; i <= nbins; i++){
		  for (int j = 1; j<= nbins;j++){
			  bias->SetBinContent(i,j,0.0);
		  }
	  }
  }
  bias->Write();


  int cameraNumber = 0;
  Int_t eventNumber;
  eventNumber = 0;

  dataset->event()->setRunNumber(run);
  (*dataset->event()->timeStamp()) = TDatime( MaxCamImageTools::getTimeStamp(infile) );

	TH2F *thisImage=0;
  int nevents=MaxCamImageTools::getNumberOfImagesInFitsFile( infile );
  cout << "Number of images = "<<nevents<<endl;
  for (int iev=0; iev<nevents; iev++) {
        dataset->event()->increaseEventNumber();


        TString tempName = "temp";
        tempName += eventNumber;
        tempName += iev; 
        thisImage = MaxCamImageTools::convertFitsIntoROOT( infile, tempName, iev);
        
        //syntax to save TH2F in a TClonesArray
        new( (*dataset->event()->ccdData())[cameraNumber]) TH2F(*thisImage);        

		int nPars=0;
		new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*driftHV);
		new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*anodeHV);
		new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*pressure);

		// ccd configuration
		new( (*dataset->event()->ccdConfig())[0]) MaxCamConfig(*ccdconfig);
        delete thisImage;
		
        dataset->fill();
		dataset->clearEventMemory();
	 
        eventNumber++;
  }

  cout << "dataset->write()" << endl;
  dataset->write();
  delete dataset;
  cout <<"Finished." << endl;
  //gSystem->Exit(0);
}


void all() {
	makeBiasFrame("data/canary_0059.root");
	for (int i=144; i<=150; i++) fitsFileProcessor(i);
}

