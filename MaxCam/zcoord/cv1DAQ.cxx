//____________________________________________________________
//
//  DAQ script for taking data.
//  Script should be run on DAQ machine only.
//  Compile as:
//
//  root> .L cv1DAQ.cxx++
//  root> init()
//  root> event(n_events)  
//  root> end()
//____________________________________________________________

#include "../DmtpcDAQ.hh"
#include "../MaxCamCamera.hh"
#include "../MaxCamWaveformTools.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>
using std::cout;
using std::endl;


DmtpcDAQ *exper;
void plotData();
int totalSaved=0; // total events
TCanvas *cwf, *cmca; 
TH1F *mca;
const int NWF=6; // max. number of waveforms to display
TString ccdType="dummy"; 
TString scopeType="ats860";
Int_t scopeNum = 0;  // identifier for the scope (only one scope in use now)
TString filename;

void init(int debug=0, TString fname="test.root", Int_t exposureTime_ms=1000);
int event(int n=1);
void end();
void go(TString fname="test.root", Int_t nevents=10, Int_t exposureTime_ms=1000);
void sequence(Int_t nruns=3, Int_t nevents=10, Int_t expTime_ms=1000, Int_t waitTime_sec=1);

TString dataDir = "/data/zcoord/";

//________________________________________
//
//    TAKE DATA
//________________________________________

void sequence(Int_t nruns, Int_t nevents, Int_t expTime_ms, Int_t waitTime_sec) {
  
  TDatime time;
  TString timestr = time.AsSQLString();
  timestr.ReplaceAll(" ", "T");
  cout << timestr << endl;

  for (Int_t ii=0; ii<nruns; ii++) {
    // make a filename
    TString fname = "sequence_";
    fname += timestr;
    fname +="_";
    fname += ii;
    fname += ".root";
    cout << fname << endl;
    // do a run
    go(fname, nevents, expTime_ms);
    // wait for some time 
    gSystem->Sleep(waitTime_sec*1000);
    // repeat
  }

}


void go(TString fname, Int_t nevents, Int_t exposureTime_ms) {
  
  TString fullFileName = dataDir+fname;
  init(0, fullFileName, exposureTime_ms);
  event(nevents);
  end();
}


void init(int debug, TString fname, Int_t exposureTime_ms) {
  filename = fname;

  exper= new DmtpcDAQ(debug,fname, ccdType, scopeType);
  int nc=exper->nCCD();

  if (ccdType!="dummy") for (int i=0; i<nc; i++) exper->ccd(i)->setTemperature(-20); 

  exper->beginRunCV1();

  for (int i=0; i<nc && ccdType!="dummy"; i++) {
    exper->ccd(i)->setVBin(4); exper->ccd(i)->setHBin(4); 

    //exper->ccd(i)->makeBiasFrame(100);
    exper->ccd(i)->readBiasFrame();
    //exper->ccd(i)->findHotPixels(100, 4, 0.1);
    //exper->ccd(i)->findHotPixels("hotpixels.dat");

    exper->ccd(i)->setNormalFrame();
  }
  
  exper->setGlobalExposureTime(exposureTime_ms);

  cwf = new TCanvas("cwf", "scope", 1200, 200);
  cwf->Divide( NWF, 1);

  cmca = new TCanvas("cmca", "scopeMCA", 400, 400);
  TString mcaTitle(fname);
  mcaTitle +=";Pulse Height (mV);#";
  mca  = new TH1F("mca", mcaTitle.Data(), 500, -250, 250);
  // make the tree readable as data rolls in
  //exper->setAutoSave(true);

  // print out the scope settings
  //exper->scope()->printScopeSetup(scopeNum);
  //cout << "input imped: " 
  //     << exper->data()->event()->scopeDataInfo(scopeNum)->getInputImpedance() << endl;

}

int event(int n) {    

  for (int i=0; i<n; i++) {
    exper->beforeEvent();
    exper->event();
    
    const char* e0="\033[44;37m";
    const char* en = "\033[0m";    
    cout << e0 << "***** EVENT SAVED ****** (n=" << totalSaved++ << ")" << en<<endl;

    // if triggered, get the pulse heights and stuff into mca
    cout << "-------- ntriggers = " << exper->data()->event()->scopeData()->GetEntries()/2 << endl;
    cout << "[" <<exper->data()->event()->scopeData()->GetEntries() << "]" <<endl;
    if (exper->data()->event()->scopeData()->GetEntries() > 0) {
      Int_t nTriggers = exper->data()->event()->scopeData()->GetEntries()/2;
      for (Int_t ii=0; ii<nTriggers; ii++) {
	MaxCamWaveformTools wfA(exper->data()->event()->scopeData(ii));
	mca->Fill(1e3*wfA.getTroughDepth());
      }
    }

    plotData();

    // clear memory after ploting
    exper->afterEvent();
  }

  return 0;
}
  
void end() {
  // save the mca to a root file
  TString outfile_mca(filename);
  outfile_mca.ReplaceAll(".root", "_mca.root");
  TFile mcaFile(outfile_mca, "RECREATE");
  mca->Write();
  mcaFile.Close();
  exper->endRun();
}

//void saveRun(char *key, char *desc) {
//  cout << key << "  " << desc << endl;
//  exper->data()->saveRun(key, desc, "NW13-039");
//}

void plotData() {

  // plot the mca
  cmca->cd(0);
  mca->Draw();
  cmca->Update();

  // plot all waveforms
  for (int iwf=0; iwf<NWF; iwf++) {
    cwf->cd(iwf+1);
    gPad->Clear();
    if (iwf>=exper->data()->event()->scopeData()->GetEntries()/2) continue;
    TH1F *hwf=exper->data()->event()->scopeData(iwf);
    //hwf->SetMaximum(1.5);
    hwf->DrawCopy();
  }
  cwf->Update();
}
