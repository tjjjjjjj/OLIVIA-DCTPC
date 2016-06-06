#include <iostream>
#include <vector>
#include <unistd.h>
#include "ScopeData.hh"
#include "ScopeWaveform.hh"
#include "ScopeHandler.hh"
#include "Scope.hh"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TTree.h"
#include "TH1F.h"  // for testing only...
#include "DmtpcDataset.hh"

#include "AlazarApi.h"
#include "AlazarCmd.h"

using std::cout;
using std::endl;
using std::vector;

int main() {

  // ensure that there are boards installed
  U32 nAlazarBoards = AlazarBoardsFound();
  if (nAlazarBoards < 1) {
    cout << "error: no alazar boards found";
    assert(!"error:  no alazar boards found");
  }

  cout << "Found " << nAlazarBoards << " board(s)" << endl;


  U32 systemId = 1;
  HANDLE h0 = AlazarGetBoardBySystemID(systemId, 0);
  HANDLE h1 = AlazarGetBoardBySystemID(systemId, 1);
  HANDLE h2 = AlazarGetBoardBySystemID(systemId, 2);

  cout << "h0, h1, h2 = " << h0 << ", " << h1 << ", " << h2 << endl;

  cout << "kind h0 = " << AlazarGetBoardKind(h0) << endl;
  cout << "kind h1 = " << AlazarGetBoardKind(h1) << endl;
  cout << "kind h2 = " << AlazarGetBoardKind(h2) << endl;


  cout << AlazarBoardsInSystemByHandle(h0) << endl;
  cout << AlazarBoardsInSystemByHandle(h1) << endl;
  cout << AlazarBoardsInSystemByHandle(h2) << endl;


  AlazarSetLED(h0, 1);
  sleep(5);
  AlazarSetLED(h0, 0);



}


//int main2() {
//
//  // First, set up a scope handler
//  // For a single scope, this is an easy way to do it:
//  //ScopeHandler* scope = new ScopeHandler(&Scope::ALAZAR_ATS860);
//  // If you want to create a scope handler and then add a scope 
//  // one-by-one, then you can do:
//  ScopeHandler* scope = new ScopeHandler();
//  scope->addScope(&Scope::ALAZAR_ATS860);
//
//  // Next, establish communication with the scope
//  scope->openScope(0);
//
//  // scopeId, state type
//  scope->configureDAQState(0, Scope::DAQ_CHARGE);
//
//  //cout << "exiting early..." << endl;
//  //return 0;
//
//  //int scopeNumber = 0;
//  //cout << "pre:  trigger level = " 
//  //     << (int)scope->getTriggerLevel1(scopeNumber) << endl;
//  //scope->setTriggerLevel1(scopeNumber,200);
//  //cout << "post: trigger level = " 
//  //     << (int)scope->getTriggerLevel1(scopeNumber) << endl;
//  //int chanNumber = 0;
//  //cout << "pre:  voltage range = " 
//  //     << (int)scope->getVoltageRange(scopeNumber,chanNumber) << endl;
//  //
//  //scope->setVoltageRange(scopeNumber, chanNumber, INPUT_RANGE_PM_10_V);
//  //cout << "post: voltage range = " 
//  //     << (int)scope->getVoltageRange(scopeNumber,chanNumber) << endl;
//
//  //cout << "return 0 to quit early" << endl;
//  //return 0;
//  
//  // Establish a ScopeData set in which triggers will be stored
//  // This must come after configureDAQState so that a ScopeConfig 
//  // instance exists
//  scope->addScopeData(0, scope->getScopeConfig(0));
//
//  // Acquire triggers for a set amount of time
//  float duration_ms = 3000.0;  // duration in milliseconds
//  int scopeNumber = 0;
//  int chanNumber = 0;
//  cout << "trigger level = " 
//       << (int)scope->getTriggerLevel1(scopeNumber) << endl;
//  cout << "voltage range = " 
//       << (int)scope->getVoltageRange(scopeNumber,chanNumber) << endl;
//
//  scope->acquireTriggers(0, duration_ms);
//
//  cout << "return 0 to quit early" << endl;
//  return 0;
//
//  // generate a 1-d vector of all waveforms generated in this event
//  vector<ScopeWaveform*> waveforms;
//  waveforms = scope->getWaveforms(0);
//
//  // set up for ROOT dataset and output
//  DmtpcDataset* _rootData = new DmtpcDataset();
//  TString outputFile = "examples/jamesdata.root";
//  _rootData->createRootFile(outputFile, "recreate"); 
//  // add a branch to the tree
//  //_rootData->tree()->Branch("data", "DmtpcDataset", &_rootData, 32000, 0);
//
//  // shove the ScopeWaveforms into the ROOT tree
//  // (but is _chargeData integrated into the tree?)
//  uint MAX_N_WF = 100;  // must match DmtpcDataset::MAX_N_TRIGGERS
//  uint ii = 0;
//  int bufSize = 32000;
//
//  bool doTH1F = true;
//  bool doWaveforms = false;
//  TClonesArray* mytca = new TClonesArray("TH1F",MAX_N_WF);
//  TClonesArray &th1f  = *mytca;
//  if (doTH1F) {
//    _rootData->tree()->Branch("branchJames", "TClonesArray", 
//			      &mytca, bufSize, 0);
//  }
//  TClonesArray* mytcawf = new TClonesArray("ScopeWaveform",MAX_N_WF);
//  TClonesArray &th1fwf  = *mytcawf;
//  if (doWaveforms) {
//    _rootData->tree()->Branch("branchJamesWf", "TClonesArray", 
//			      &mytcawf, bufSize, 0);
//  }
//
//  /* save this for reference
//    // a way to loop over a TClonesArray
//    for (int i=0; i<array->GetEntriesFast(); i++) {
//      MyClass* p = (MyClass*) array->At(i);
//      p->methodOfMyClass();
//    }
//  */
//  int chAId = 0;
//  int nbins = int(scope->getScopeConfig(0)->board()->getRecordSize());
//  ScopeWaveform*  wfChanA;
//  ScopeWaveform* myswf;
//  TH1F* myhist;
//  //while ( (ii<3 ) ) {
//  while ( (ii<waveforms.size()) && (ii<MAX_N_WF) ) {
//    TString iname = "histname";
//    iname += ii;
//    // populate the TClonesArray
//    // use a copy constructor...
//    if (doWaveforms) 
//      myswf = new(th1fwf[ii]) ScopeWaveform( *(scope->data(0)->dataChan(chAId)->wf(ii)) );
//    if (doTH1F)
//      myhist = new(th1f[ii]) TH1F(iname, "hist title", nbins, 0, nbins);
//
//    wfChanA = scope->data(0)->dataChan(chAId)->wf(ii);
//
//    // get access to the object stored in the TClonesArray.
//    // notice:
//    // (1) you can use the pointer to TCA
//    // (2) you must cast the pointer returned by At().
//    // (3) you must use At(), you cannot just index into the TCA 
//    //     (e.g. mytcawf[ii] does not work).
//    // populate the histogram
//    for (int jj=0; jj<nbins; jj++) {
//      if (doWaveforms)
//	myswf->SetBinContent(jj, (ii+1)*jj);
//      if (doTH1F)
//	myhist->SetBinContent(jj, float(wfChanA->GetBinContent(jj)));
//    }
//    ii++;
//  }
//  
//  _rootData->fill();
//  _rootData->write();
//  //_rootData->tree()->Print();
//
//  scope->clearData(0);
//  th1f.Delete();
//
//  //  if (scope->wasTriggered()) {
//  //    // read out camera
//  //    // save camera data to file
//  //    // save scope waveform data to file
//  //  }
//  return 0;
//}

