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
#include "../MaxCamConfig.hh"
#include "../MaxCamImageTools.hh"
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


const int NCAM = 2;
//                       A         B         C        D
//TString sns[NCAM] = {"081013", "100534", "100439", "081264"};
//TString sns[NCAM] = {"100534", "100439", "081264"};
//TString rotangs[NCAM] = {"QUARTER_LEFT", "QUARTER_RIGHT", "HALF_TURN"};
TString sns[NCAM]     = {"100534", "100439"};
TString rotangs[NCAM] = {"QUARTER_LEFT", "QUARTER_RIGHT"};
//TString sns[NCAM] = {"081264"};
//TString rotangs[NCAM] = {"HALF_TURN"};
//TString rotangs[NCAM] = {"NONE", "QUARTER_LEFT", "QUARTER_RIGHT", "HALF_TURN"};
//TString rotangs[NCAM] = {"HALF_TURN", "QUARTER_LEFT", "QUARTER_RIGHT", "NONE"};
bool doRotate = false;

DmtpcRun *exper;
void plotCCD();
void plotScope();
void plotData();
int totalSaved=0, totalAll=0; // total events
TCanvas *c=0, *cwf; 
const int NWF=6; // max. number of waveforms to display + sum
TH1F *hmca=0, *hmcb=0;
TString ccdType="dummy"; 
long globalExpoTime=1000; //msec
long binning=4;
TString scopeType="alazar";
bool doPlot=false;
bool usePulseMaximum=true;
int adcBits=16;

bool verbose=false;

//Use these to prevent the DAQ process from working more than once at a time
bool createLockFile(const char * name)
{
        int fd = open(name, O_WRONLY | O_CREAT | O_EXCL); 
	if (fd < 0) return false; 
	close(fd); 
	return true; 
}


void releaseLockFile(const char * name)
{
  unlink(name); 
}

//________________________________________
//
//    TAKE CCD DATA
//________________________________________

void init(int debug=0, TString fname="test.root", TString detid="", Int_t nbiasframes=100 ) {

  exper= new DmtpcRun(debug,fname, ccdType, scopeType);
  int nc=exper->nCCD();

  exper->data()->setDetId(detid); 

  if (ccdType!="dummy") for (int i=0; i<nc; i++) exper->ccd(i)->setTemperature(-20); 
  exper->beginRun();


  exper->setGlobalExposureTime(globalExpoTime); 

  for (int i=0; i<nc && ccdType!="dummy"; i++) {
    cout << "setting digitize overscan." << endl;
    exper->ccd(i)->setDigitizeOverscan(true);
    cout << "setting binning" << endl;
    //cout << "exper->ccd("<<i<<")->setVBin("<<binning<<")"<<endl;
    exper->ccd(i)->setVBin(binning); 
    //cout << "exper->ccd("<<i<<")->setHBin("<<binning<<")"<<endl;
    exper->ccd(i)->setHBin(binning); 
    //cout << "exper->ccd("<<i<<")->makeBiasFrame("<<nbiasframes<<")" << endl;
    //exper->ccd(i)->makeBiasFrame(nbiasframes);
    //exper->ccd(i)->readBiasFrame();
    cout <<"Setting ADC bits to "<<adcBits<<endl;
    exper->ccd(i)->setDataBits(adcBits);
  }
  exper->makeBiasFrames(nbiasframes);
  
  gStyle->SetOptTitle(kTRUE);

}

void initPlot() {
  if (verbose) cout << "initPlot() -----------------" << endl;
  c = new TCanvas("c", "CCD", 0,0, 810, 550);


  if (NCAM == 2) {
    c->Divide( 2, 1); // need to generalize this to N cameras...
  } else {
    c->Divide( 2, 2); // need to generalize this to N cameras...
  }

  cwf = new TCanvas("cwf", "Scope", 820, 0, 440, 900);
  int nscope = exper->nScope();
  cwf->Divide( 2*nscope, NWF);
  cout << "the waveform canvas has N subpads.  N = " << cwf->GetListOfPrimitives()->GetEntries() << endl;
  hmca =new TH1F("hmca", "",128, 0-0.03, exper->getScopeVoltageMax(0,0)-0.03);
  hmcb =new TH1F("hmcb","", 128, 0-0.03, exper->getScopeVoltageMax(0,1)-0.03);
}


void setDriftVoltage(double val) {
  assert(val>-1e-5);
  exper->setDriftVoltage(val);
}

void setAnodeVoltage(double val) {
  assert(val>-1e-5&&val<1.5);
  exper->setAnodeVoltage(val);
}

void setPressureAndSave(double val) {
  assert(val>0);
  exper->setPressureAndSave(val);
}



//
//float yieldMin=-50, yieldMax=250;
float yieldMin=-50, yieldMax=1000;
int nrecovery=0;
int event(int n=1) {    
  
  if (doPlot && c==0) initPlot();
  
  for (int i=0; i<n; i++) {
    
    totalAll++;
    exper->beforeEvent();
    
    cout << "dmtpcDAQ::event() pre exper->event()" << endl;
    exper->event();
    cout << "dmtpcDAQ::event() post exper->event()" << endl;
    
    if (exper->getSaveFlag()) {
      const char* e0="\033[44;37m";
      const char* en = "\033[0m";    
      cout <<e0 << "***** EVENT SAVED ****** (n=" << totalSaved++ << "/" << n << ")" << en<<endl;
    }
    
    cout << "doPlot = " << doPlot << endl;
    if (doPlot) plotData();
    
    // clear memory after ploting
    cout << "dmtpcDAQ::event() pre exper->afterEvent()" << endl;
    exper->afterEvent();
  }
  
  
  return 0;
}


int focusEvent(int nev=100) {
  for (int iccd=0; iccd<exper->nCCD(); iccd++) { exper->ccd(iccd)->openShutter(); }
  for (int i=0; i<nev; i++) { event(); c->Update(); }
  for (int iccd=0; iccd<exper->nCCD(); iccd++) { exper->ccd(iccd)->closeShutter(); }
}



void end() {
  exper->endRun();
}


void saveRun(char *scpdest="mitdm00.mit.edu:/data/") {
  exper->data()->saveRun(scpdest);
}


void plotData() {
  plotCCD();
  plotScope();
}

void plotScope() {
  int nscope = exper->nScope();
  int nchan = nscope*2;

  int nwf = exper->data()->event()->scopeData()->GetEntries();
  int ntrg = nwf/nchan;

  cout << "nscope, nwf, ntrg = " << nscope << ", " << nwf << ", " << ntrg << endl;
  // triggers are stored as [ch1_trg1, ... ,ch1_trgN, ch2_trg1, ..., ch2_trgN, ...]
  cwf->Clear("D");  // clear plots but preserve sub-pads
  int ipad,index;
  for (int itrg=0; itrg<ntrg; itrg++) {
    cout << "itrg = " << itrg << endl;
    for (int ich=0; ich<nchan; ich++) {
      ipad = (itrg*nchan)+ich+1;
      index = itrg + ich*ntrg;
      cout << "  ipad, index = " << ipad << ", " << index << endl;
      cout << "      exper->data()->event()->scopeData()->GetEntries() = " << exper->data()->event()->scopeData()->GetEntries() << endl;
      cwf->cd(ipad);
      exper->data()->event()->scopeData(index)->DrawCopy();
    }
  }

  cwf->Update();
}

void plotCCD() {
  // plot all CCD images

  int ipad=0;
  if (verbose) cout << "plotData()" << endl;
  for (int iccd=0; iccd<exper->data()->event()->ccdData()->GetEntries()&&ccdType!="dummy"; iccd++) {

    cout << endl;
    cout << endl;
    cout << "iccd = " << iccd << endl;
    
    for (int ii=0; ii<NCAM; ii++) {
      if (sns[ii] == exper->data()->event()->ccdConfig(iccd)->serialNumber) {
	ipad = ii+1;
	cout << "ii, sns[ii], ipad = " << ii << ", " << sns[ii] << ", " << ipad << endl;
	break;
      }
    }
    
    c->cd(ipad);

    TH2F *hccd=0;
    if (exper->data()->event()->ccdData()->GetEntries()) {
      hccd=exper->data()->event()->ccdData(iccd);
      cout << "exper->data()->event()->ccdData("<<iccd<<") = " << exper->data()->event()->ccdData(iccd) << endl;
      cout << "exper->data()->event()->ccdData("<<iccd<<")->GetEntries() = " << exper->data()->event()->ccdData(iccd)->GetEntries() << endl;
      
      /////////////////////////////////////
      // start kluge
      // for some reason, histos do not plot with "colz" option... they *do* with "lego" option
      // GetEntries() on the histogram returns zero!
      // so forcibly repopulate the histogram with its values!?!? it works...
      cout << "WARNING:  KLUGE IN PLACE FOR PLOTTING IMAGES..." << endl;
      cout << "WARNING:  BINNING = 4 IS ASSUMED" << endl;
      for (int ii=0; ii<257; ii++) {
      	for (int jj=0; jj<257; jj++) {
      	  hccd->SetBinContent(ii, jj, hccd->GetBinContent(ii, jj));
      	}
      }
      cout << "exper->data()->event()->ccdConfig("<<iccd<<")->serialNumber = " << exper->data()->event()->ccdConfig(iccd)->serialNumber << endl;
      hccd->SetTitle(exper->data()->event()->ccdConfig(iccd)->serialNumber);
      // end kluge
      /////////////////////////////////////

      if (exper->ccd(iccd)->biasFrame()&&hccd) {
	cout << "found bias frame ... subtracting" << endl;
	hccd->Add( exper->ccd(iccd)->biasFrame(), -1);
      }
      hccd->SetMaximum(yieldMax);
      hccd->SetMinimum(yieldMin);
      if (doRotate) {
	cout << "rotatePerfect" << endl;
	MaxCamImageTools::rotatePerfect(hccd, rotangs[ipad-1])->DrawCopy("colz");
      } else {
	hccd->DrawCopy("colz");
      }
    }
    
    c->Update();
  }

}
