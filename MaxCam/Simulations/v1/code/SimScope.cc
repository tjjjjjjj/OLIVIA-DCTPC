// MaxCam includes
#include "SimScope.hh"

// ROOT includes
#include "TObject.h"

// C++ includes
#include <iostream>
using std::cout;
using std::endl;
#include <assert.h>

//________________________________
/*Begin_Html
<center><h2>The SimScope Class </h2></center>
The SimScope class is designed to hold the various parameters of the 
Scopes used in the detector.  A separate class is used so that in
future versions of the simulation, several instances of the class may
be added to a TObjectArray or TClonesArray to allow for generation of
multiple Scope signal, each having its own properties.
End_Html*/
//________________________________

ClassImp(SimScope)

SimScope::SimScope()
{
  //<<<<<Default Constructor>>>>>
  fSerialNumber = "";
  fScopeNumber = 0;
  fChans = new TObjArray();
  fWfs = new TObjArray();
}
SimScope::SimScope(TString ser, Int_t pmtnum)
{
  //<<<<<Alternate Constructor>>>>
  fSerialNumber = ser;
  fScopeNumber = pmtnum;
  fChans = new TObjArray();
}

SimScope::~SimScope()
{
  //<<<<<Destructor>>>>>
}

void
SimScope::setBoardType(TString bt) {
  fBoardType= bt;
  updateHardwareSpecs();
}

void
SimScope::updateHardwareSpecs() {
  // given the board type, set all of the hardware parameters
  if (fBoardType=="ALAZAR_ATS860") {
    cout << "found an Alazar ATS860 board" << endl;
    fNChan = 2;
    fNbitsPerSamp = 8;
  }
}

Float_t
SimScope::getVoltageStep(Int_t ch)
{
  Int_t nlevels = getNlevels();
  Float_t vmin  = chan(ch)->getVoltageMin();
  Float_t vmax  = chan(ch)->getVoltageMax();
  return (vmax-vmin)/nlevels;
}

Float_t
SimScope::getLevelOfVolts(Int_t ch, Float_t level_volts) 
{
  Float_t vmin  = chan(ch)->getVoltageMin();
  int level     = (level_volts-vmin)/getVoltageStep(ch);

  Int_t nlevels = getNlevels();
  if (level<0) level=0;
  if (level > nlevels) level=nlevels-1;
    
  return level;
}


ScopeDataInfo* 
SimScope::chan(Int_t channum) {
  int nchan = fChans->GetEntries(); 
  ScopeDataInfo* chan;
  if (channum < nchan) chan = (ScopeDataInfo*) fChans->At(channum);
  else{
    cout << "Invalid channel number: "<<channum<<". There are "<<nchan <<" channels. Returning null pointer." <<endl;
    chan = NULL;
  }
  return chan;
  
}



//Used for debugging.
int checks = 0;
void
check(int& checks)
{
  cout << "Checkpoint " << checks << endl;
  checks++;
}

void
SimScope::setWfs(SimScope *scope, Int_t channum) {
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
  cout << "scope wf name = " << name << endl;
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
SimScope::resetWfs() {
  for (Int_t ii=0; ii<fWfs->GetEntries(); ii++) 
    ((ScopeWaveformData*)fWfs->At(ii))->Reset();
}

ScopeWaveformData*
SimScope::wf(Int_t itrig) 
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
