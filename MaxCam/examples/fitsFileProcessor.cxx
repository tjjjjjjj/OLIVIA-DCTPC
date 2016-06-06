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

// function prototypes
void fitsFileProcessor(TString fitslist="files.dat", TString outputFile="dmtpc_run00000.root", TString darkFitsFile="dark.fits");
TObjArray* readFitsFilelist(TString listname="files.dat");


void fitsFileProcessor(TString fitslist, TString outputFile, TString darkFitsFile) {

  cout << "fitslist=" << fitslist << endl;
  cout << "outputFile=" << outputFile << endl;
  
  // MOST OF THESE NUMBERS ARE JUST DUMMY NUMBERS
  
  // inputs
  //int npix=1024;
  int npix=512;
  int nbins=512;
  //float noise=7;

  cout << "inside fitsFileProcessor" << endl;

  cout << "Initializing..."<<endl;
  
  ///////////////////////////////////////////////////////////////////////////////
  ////Set up data file                                                         //
  ///////////////////////////////////////////////////////////////////////////////
  DmtpcDataset *dataset = new DmtpcDataset();
  dataset->createRootFile(outputFile,"recreate");
  
  dataset->setComment("COMMENT");
  dataset->setLocation("LOCATION");
  dataset->setKeyword("KEYWORD");
  dataset->comment()->Write();
  dataset->keyword()->Write();
  dataset->location()->Write();

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

  //////////////////////////////////////////////////////////////////////////////////
  ////Set up MaxCamConfig file for the dataset.                                   //
  //////////////////////////////////////////////////////////////////////////////////
  cout << "Set up MaxCamConfig" << endl;
  
  MaxCamConfig *ccdconfig = new MaxCamConfig("ccdConfig","CCD Camera configuration");
  ccdconfig->cameraID = 0;
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
  
  // populating the dark frame
  //TH2F *bias = new TH2F("biasFrame1","biasFrame1",nbins,0,nbins,nbins,0,nbins);
  //for (int i = 1; i <= nbins; i++){
  //  for (int j = 1; j<= nbins;j++){
  //    //bias->SetBinContent(i,j,gRandom->Gaus(0,noise));
  //    // empty bias frame
  //    bias->SetBinContent(i,j,0.0);
  //  }
  //}
  TH2F* bias = MaxCamImageTools::convertFitsIntoROOT(darkFitsFile,"biasFrame1", 0);
  bias->Write();

  //////////////////////////////////////////////////////////////
  ////Save CCD image and fill dataset                         //
  //////////////////////////////////////////////////////////////
  cout << fitslist << endl;
  TObjArray* fitsFiles = readFitsFilelist(fitslist);
  TIter next(fitsFiles);
  while (TObjString* fitsFile=(TObjString *)next()) {
    int cameraNumber = 0;
    Int_t eventNumber;
    eventNumber = 0;
    cout << "[" << fitsFile->GetString() << "]" << endl;
  
    // see DmtpcDAQ::beforeEvent()
    dataset->event()->setRunNumber(1);
    // ignore the timestamp for now (the ProEM fits files don't have a DATE field in the FITS header)
    //(*dataset->event()->timeStamp()) = TDatime(MaxCamImageTools::getTimeStamp(fitsFile->GetString()));
    //Events start from 1.
    dataset->event()->increaseEventNumber();
  
    int nevents=MaxCamImageTools::getNumberOfImagesInFitsFile(fitsFile->GetString());
    //cout << "Number of images = "<<nevents<<endl;
    for (int iev=0; iev<nevents; iev++) {
      TString tempName = "temp";
      tempName += eventNumber;
      tempName += iev; 
      //cout << tempName << endl;
      TH2F* thisImage = MaxCamImageTools::convertFitsIntoROOT(fitsFile->GetString(),tempName, iev);
        
      //syntax to save TH2F in a TClonesArray
      new( (*dataset->event()->ccdData())[cameraNumber]) TH2F(*thisImage);        

      int nPars=0;
      new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*driftHV);
      new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*anodeHV);
      new( (*dataset->event()->experimentConfig())[nPars++] ) MaxCamChannel(*pressure);

      // ccd configuration
      new( (*dataset->event()->ccdConfig())[0]) MaxCamConfig(*ccdconfig);

      // must call "delete" *after* dataset->fill()
      dataset->fill();
      delete thisImage;
      dataset->clearEventMemory();
      eventNumber++;
    }

  }

  cout << "dataset->write()" << endl;
  dataset->write();
  cout << "dataset->file()->Close()" << endl;
  dataset->file()->Close();
  //delete dataset;
  cout <<"Finished." << endl;
  //gSystem->Exit(0);
}



TObjArray* readFitsFilelist(TString listname) {
  // read a list of fits files from an ASCII file

  // to hold the list of fits files to be 
  // converted into a single DmtpcDataset root file
  TObjArray *fitsList = new TObjArray();

  Bool_t verbose = 0;
  // read in the list of fits files
  ifstream ifs(listname.Data());
  string line;
  while (getline(ifs,line)) {
    if (verbose) cout << "[" << line << "]" << endl;
    TObjString *filename = new TObjString(line.c_str());
    fitsList->Add(filename);
  }

  if (verbose) {
    TIter next(fitsList);
    while (TObjString* myFitsFile=(TObjString *)next()) {
      //cout << fitsList->At(ii) << endl;
      cout << myFitsFile->GetString() << endl;
    }
  }

  return fitsList;
}
