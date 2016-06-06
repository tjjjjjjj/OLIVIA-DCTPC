// this is built off of looking at Jeremy's makeSimulations.cc

// headers
#include "TCanvas.h"
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
#include "../SPEFile.hh"

#include "../DmtpcDataset.hh"
#include "../MaxCamConfig.hh"
#include "../MaxCamImageTools.hh"
#include "../MaxCamChannel.hh"

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

// function prototypes
// In SPEFileProcessor(): if outputlist=="", each filename.SPE file is automatically saved to filename.root file in same directory. Otherwise, specify name of list of SPE output files in outputlist.
void SPEFileProcessor(TString spelist="files.dat", TString outputlist="", TString darkFitsFile="dark.fits");
TObjArray* readSPEFilelist(TString listname="files.dat");

void SPEFileProcessor(TString spelist, TString outputlist, TString darkFitsFile) {

  cout << "spelist=" << spelist << endl;
  cout << "outputlist=" << outputlist << endl;
  
  // MOST OF THESE NUMBERS ARE JUST DUMMY NUMBERS
  
  // inputs
  //int npix=1024;
  int npix=512;
  int nbins=512;
  //float noise=7;

  cout << "inside SPEFileProcessor" << endl;

  cout << "Initializing..."<<endl;
  


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
  
  //// populating the dark frame
  //TH2F *bias = new TH2F("biasFrame1","biasFrame1",nbins,0,nbins,nbins,0,nbins);
  //for (int i = 1; i <= nbins; i++){
  //  for (int j = 1; j<= nbins;j++){
  //    //bias->SetBinContent(i,j,gRandom->Gaus(0,noise));
  //    // empty bias frame
  //    bias->SetBinContent(i,j,0.0);
  //  }
  //}
  ////TH2F* bias = MaxCamImageTools::convertFitsIntoROOT(darkFitsFile,"biasFrame1", 0);
  //bias->Write();

  //////////////////////////////////////////////////////////////
  ////Save CCD image and fill dataset                         //
  //////////////////////////////////////////////////////////////
  cout << spelist << endl;
  TObjArray* speFiles = readSPEFilelist(spelist);
  TIter next(speFiles);
  DmtpcDataset *dataset = new DmtpcDataset();

  TString rootOutDir = "/export/data03/jbattat/cameras/proem/root/";
  
  while (TObjString* speFilen=(TObjString *)next()) {
    int cameraNumber = 0;
    Int_t eventNumber;
    eventNumber = 0;
    TString filename=speFilen->GetString();

    //Set up data file   
    //if (outputlist==""){
    //}
    //else {
    //  TObjArray* outFiles = readSPEFilelist(outputlist);
    //  TIter next(outFiles);
    //  TObjString* outFilen=(TObjString *)next();
    //  dataset->createRootFile(outFilen->GetString(), "recreate");
    //}
    
    int strlen=filename.Length();
    int slashindex=filename.Last('/')+1;
    TString rootFile = filename(slashindex, strlen);
    rootFile.ReplaceAll(".SPE", ".root");
    TString rootFileOut = rootOutDir + rootFile;

    cout << "speFileIn   = " << filename << endl;
    cout << "rootFileOut = " << rootFileOut << endl;

    dataset->createRootFile(rootFileOut,"recreate");

    dataset->setComment("COMMENT");
    dataset->setLocation("LOCATION");
    dataset->setKeyword("KEYWORD");
    dataset->comment()->Write();
    dataset->keyword()->Write();
    dataset->location()->Write();

    // see DmtpcDAQ::beforeEvent()
    dataset->event()->setRunNumber(1);

    // ignore the timestamp for now (the ProEM fits files don't have a DATE field in the FITS header)
    //(*dataset->event()->timeStamp()) = TDatime(MaxCamImageTools::getTimeStamp(speFile->GetString()));

    SPEFile imagefile(filename);

    int nevents=imagefile.Getnimages();
    cout << "Number of images = "<<nevents<<endl;
    for (int iev=0; iev<nevents; iev++) {
    //TCanvas *c1 = new TCanvas();
      //for (int iev=0; iev<1; iev++) {
    //Events start from 1.
    dataset->event()->increaseEventNumber();
    TH2F* thisImage = imagefile.GetImagetoROOT(iev);
    //thisImage->Draw("colz");
    //c1->Update();
  
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
    } // loop over images
    
    cout << "dataset->write()" << endl;
    dataset->write();
    cout << "dataset->file()->Close()" << endl;
    dataset->file()->Close();
    //delete dataset;
    cout <<"Finished." << endl;
    //gSystem->Exit(0);
  } // loop over files
  
}




TObjArray* readSPEFilelist(TString listname) {
  // read a list of spe files from an ASCII file

  // to hold the list of spe files to be 
  // converted into a single DmtpcDataset root file
  TObjArray *speList = new TObjArray();

  Bool_t verbose = 0;
  // read in the list of spe files
  ifstream ifs(listname.Data());
  string line;
  while (getline(ifs,line)) {
    if (verbose) cout << "[" << line << "]" << endl;
    TObjString *filename = new TObjString(line.c_str());
    speList->Add(filename);
  }

  if (verbose) {
    TIter next(speList);
    while (TObjString* myspeFile=(TObjString *)next()) {
      //cout << speList->At(ii) << endl;
      cout << myspeFile->GetString() << endl;
    }
  }

  return speList;
}





