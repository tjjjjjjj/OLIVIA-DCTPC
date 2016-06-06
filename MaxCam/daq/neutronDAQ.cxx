//____________________________________________________________
//
//  DAQ script for taking data.
//  Script should be run on DAQ machine only.
//  Compile as:
//
//  root> .L dmtpcDAQ.cxx++
//  root> init()
//  root> event(n_events)  
//____________________________________________________________

#include "../DmtpcRun.hh"
#include "../MaxCamCamera.hh"
#include "../ScopeDataInfo.hh"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include <fcntl.h>
#include <stdio.h>
#include "TStyle.h"

#include <iostream>
using std::cout;
using std::endl;


DmtpcRun *exper;
void plotData();
int totalSaved=0, totalAll=0; // total events
TCanvas *c=0, *cwf; 
const int NWF=6; // max. number of waveforms to display + sum
TH1F *hmca=0, *hmcb=0;
TString ccdType="dummy"; 
long globalExpoTime=1000; //msec
long binning=1;
TString scopeType="alazar";
bool doPlot=true;
bool usePulseMaximum=true;
float ccdTemp=-20;
int nBiasFrames=100; 
int ccdGain=0;


//________________________________________
//
//    TAKE CCD DATA
//________________________________________


void init(int debug=0, TString fname="test.root", TString detid="Raytheon") {


  exper= new DmtpcRun(debug,fname, ccdType, scopeType);
  int nc=exper->nCCD();
  exper->data()->setDetId(detid); 

  if (ccdType!="dummy") for (int i=0; i<nc; i++) exper->ccd(i)->setTemperature(ccdTemp); 
  exper->beginRun();


  exper->setGlobalExposureTime(globalExpoTime);
  for (int i=0; i<nc && ccdType!="dummy"; i++) {
    exper->ccd(i)->setVBin(binning); exper->ccd(i)->setHBin(binning); 
    exper->ccd(i)->setGain(ccdGain);
    exper->ccd(i)->makeBiasFrame(nBiasFrames);
    //exper->ccd(i)->readBiasFrame();
  }
  
  gStyle->SetOptTitle(kTRUE);

}

void initPlot() {
 
  c = new TCanvas("c", "CCD", 0,0, 600, 500);
  cwf = new TCanvas("cwf", "CCD", 820, 0, 440, 900);
  int nscope = exper->nScope()>0 ? exper->nScope() : 1;
  cwf->Divide( 2*nscope, NWF);

  //cout << "Init plot done" << endl;
}


void setDriftVoltage(double val) {
  exper->setDriftVoltage(val, "drift_hv");
}

void setAnodeVoltage(double val) {
  exper->setAnodeVoltage(val, "anode_hv");
}

void setPressureAndSave(double val) {
  assert(val>0);
  exper->setPressureAndSave(val);
}

void rampDown() {

}
void rampUp() {

}


//
float yieldMin=-5, yieldMax=5;
int nrecovery=0;
int event(int n=1) {    
  
  if (doPlot && c==0) initPlot();

  for (int i=0; i<n; i++) {
    
    totalAll++;
    exper->beforeEvent();

    exper->event();

    if (exper->getSaveFlag()) {
      const char* e0="\033[44;37m";
      const char* en = "\033[0m";    
      cout <<e0 << "***** EVENT SAVED ****** (n=" << ++totalSaved << "/" << totalAll << ")" << en<<endl;
    }

    if (doPlot) plotData();
    
    // clear memory after ploting
    exper->afterEvent();
  }


  return 0;
}
  

void end() {
  exper->endRun();
}


void saveRun() {
  exper->data()->saveRun("/data/");
}




void plotData() {
  // plot all CCD images
  for (int iccd=0; iccd<exper->data()->event()->ccdData()->GetEntries()&&ccdType!="dummy"; iccd++) {
    c->cd();
    if (exper->data()->event()->ccdData()->GetEntries()) {
      TH2S *hccd=(TH2S*)exper->data()->event()->rawCcdData(iccd);
      if ( exper->ccd(iccd)->biasFrame() && hccd ) hccd->Add( exper->ccd(iccd)->biasFrame(), -1);
      hccd->SetMaximum(yieldMax);
      hccd->SetMinimum(yieldMin);
      hccd->DrawCopy("colz");
    }
  }
  // plot waveforms
  int nscope=exper->nScope();
  for (int i=1; i<=NWF*2*nscope; i++) {
    cwf->cd(i); gPad->Clear();
  }
  int nwf=exper->data()->event()->scopeData()->GetEntries();
  for (int iwf=0; iwf<nwf; iwf++) {
    TH1F *hwf=exper->data()->event()->scopeData(iwf);
    const char *hname = hwf->GetName();
    
    assert(hname);
    int boardID, event;
    char channel;
    sscanf(hname, "scope_%d_%c_%d", &boardID, &channel, &event);
    int channelID = channel=='A' ? 0 : 1;

    int ipad = nscope*2*event + boardID*2 + channelID +1;
    cwf->cd(ipad);
    hwf->DrawCopy();
    

    //cout << hname << ": "<< boardID << ",  " << channelID << ",  " << event << endl;

  }

  //cout << "total wf="<<nwf << " plotted="<<iwf<<endl;
}
