#include "SimPMT.hh"

#include "../../../ScopeWaveformData.hh"

//#include "TH1.h"
#include "TObject.h"
//#include "TF1.h" 
//#include "TFile.h"
//#include "TRandom3.h"
#include <iostream>
#include <math.h>
#include <assert.h>
using std::cout;
using std::endl;
//________________________________
/*Begin_Html
<center><h2>The SimPMT Class </h2></center>
The SimPMT class is designed to hold the various parameters of the 
PMTs used in the detector.  A separate class is used so that in
future versions of the simulation, several instances of the class may
be added to a TObjectArray or TClonesArray to allow for generation of
multiple PMT signal, each having its own properties.
End_Html*/
//________________________________

ClassImp(SimPMT)

SimPMT::SimPMT()
{
  //<<<<<Default Constructor>>>>>
  fSerialNumber = "";
  fPMTNumber = 0;
  fX = 0;
  fY = 0;
  fDiameter = 25.4;
  fGain  = 1e4;
  fWlMin  = 200.;
  fWlMax  = 1000.;
  fWlStep = 1.;
  fPMTQE = NULL;
  setQE("flat",1.0);
  //setFlatQE(1.0);
  fWfs = new TObjArray();
  fScope=NULL;
  fScopeChan=-1;
}

SimPMT::SimPMT(double x, double y,  double diam, double gain, int pmtnum)
{
  //<<<<<Alternate Constructor>>>>
  fSerialNumber = "";
  fPMTNumber = pmtnum;
  fX = x;
  fY = y;
  fDiameter = diam;
  fGain = gain;
  fWlMin  = 200.;
  fWlMax  = 1000.;
  fWlStep = 1.;
  fPMTQE = NULL;
  setQE("flat",1.0);
  //setFlatQE(1.0);
  fWfs = new TObjArray();
  fScope=NULL;
  fScopeChan=-1;
}

SimPMT::~SimPMT()
{
  //<<<<<Destructor>>>>>
}

ScopeWaveformData*
SimPMT::wf(Int_t itrig) 
{
  int nwf = fWfs->GetEntries();
  ScopeWaveformData* wf;
  if (itrig < nwf) wf = (ScopeWaveformData*) fWfs->At(itrig);
  else {
    cout << "Invalid WF number: "<<itrig<<".  There are "<<nwf<<" triggers.  Returning null pointer." << endl;
   wf = NULL;
  }
  return wf;
}

void
SimPMT::resetWfs() {
  for (Int_t ii=0; ii<fWfs->GetEntries(); ii++) 
    ((ScopeWaveformData*)fWfs->At(ii))->Reset();
}

void
SimPMT::setWfs(SimScope *scope, Int_t channum) {
  //Sets the PMT waveform according to the associated SimScope channel
  //configuration

  Int_t scopenum = scope->getScopeNumber();

  /* FIXME: hardcoded triggerTime */
  Int_t triggerTime = 0;

  /* FIXME: hardcoded channel letter in histname, should be determined by channum 0=A, 1=B, etc. */
  TString chanletter;
  if (channum==0) chanletter="A";
  else if (channum==1) chanletter="B";
  else {
    cout << "only prepared for a two channel scope...  " << endl;
    cout << "you gave channum = " << channum << endl;
    assert(!"invalid channel number");
  }
  // See DmtpcEvent::scopeData(trigger, board, channel) for wf naming convention
  /* FIXME:  hardcoded trigger number...*/
  Int_t trignum=1; // for some annoying reason, trignum is 1-indexed, not 0-indexed.
  TString name = TString::Format("scope_%d_%s_%d",scopenum, chanletter.Data(),trignum);
  cout << "pmt wf name = " << name << endl;

  Int_t nbins, nbinsPre;
  Float_t xlow, xup, levelsToVolts, dtScope, vZero;

  levelsToVolts = scope->getVoltageStep(channum);
  vZero         = float(scope->getLevelOfVolts(channum, 0.0));
  nbins         = scope->getRecordLength();
  nbinsPre      = scope->getRecordPreSize();
  dtScope       = 1.0/scope->getClockRate();
  xlow = -nbinsPre*dtScope;
  xup  = (nbins-nbinsPre)*dtScope;

  /* FIXME:  add check to see if fWfs is empty... */
  // the idea here is to initialize fWfs with one ScopeWaveformData entry 
  // if there are multiple triggers per event, then this will need to 
  // be extended...
  fWfs->Add(new ScopeWaveformData(name, nbins, xlow, xup, triggerTime, levelsToVolts, vZero));

  for (int i=0; i<fWfs->GetEntries(); i++)
    ((ScopeWaveformData*)fWfs->At(i))->GetXaxis()->SetTitle("time [ns]");
}


void
SimPMT::setPosition(double x, double y)
{
  //Sets the x-y position of the image center
  fX = x;
  fY = y;
}

void
SimPMT::setQE(TString filename, double val, TRandom3 *rnd) {
  //if (fPMTQE == NULL) 
    fPMTQE = new DmtpcPMTQE(filename, val, rnd);
}

