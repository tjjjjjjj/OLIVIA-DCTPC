#include <iostream>
#include <string>
#include <vector>

#include "TClonesArray.h"

#include "ScopeHandler.hh"
#include "Scope.hh"
#include "ScopeAlazarATS860.hh"
#include "ScopeRecord.hh"

using std::cout;
using std::endl;
using std::string;
using std::vector;

int main() {

  cout << "main()\n";
  
  ScopeHandler *sh;
  sh = new ScopeHandler(Scope::ALAZAR_ATS860);

  // can access the ith scope with sh->scope(i)
  // or sh->scope() gives 0th scope

  cout << "sh->scope()->openScope(0)\n";
  sh->scope()->openScope(0);

  cout << "sh->makeScopeConfig(ScopeHandler::SCOPE_SETUP_TEST)" << endl;
  sh->scope()->makeScopeConfig(Scope::SCOPE_SETUP_TEST);

  cout << "sh->scope()->initScope()\n";
  if (sh->scope()->initScope()) {
    cout << "failed scope initialization" << endl;
    return -1;
  }

  int time0 = time(0);
  cout << "sh->scope()->acquireTriggers()\n";
  sh->scope()->acquireTriggers(1);
  int time1 = time(0);
  cout << "  -- actually took " << (time1-time0) << " seconds" << endl;

  cout << "sh->scope()->getNValidTriggers()\n";
  cout << "There were [" << sh->scope()->getNValidTriggers()
       << "] valid triggers." << endl;

  vector<ScopeRecord*>* myRecords;
  myRecords = sh->scope()->getRecords(1);

  //sh->scope()->createHistosFromWaveforms();
  int nTriggers = sh->scope()->getNValidTriggers();
  if (nTriggers <= 0) return 0;


  int nmax = 100;
  TClonesArray* chargeData;
  //chargeData = new TClonesArray("TScopeRecord", nmax);
  chargeData = new TClonesArray("TH1F", nmax);
  int iSavedWf=0;
  while ( (iSavedWf < nTriggers) && (iSavedWf < nmax) ) {
    //(*chargeData)[iSavedWf] = sh->scope()->recordOfWaveform(iSavedWf);
    //(*chargeData)[iSavedWf] = new TH1F("name", "", 100, 0, 100);
    (*chargeData)[iSavedWf] = sh->scope()->histoOfWaveform(iSavedWf);
    iSavedWf++;
  }

  return 0;

  
  (*myRecords->begin())->getSampleRate();
  //cout << "sample rate = " << (int)myRecords->first()->getSampleRate(); << endl;



  cout << "sh->scope()->getChannelData(1)" << endl;
  vector<ScopeWaveform*>* tester;
  tester = sh->scope()->getChannelData(1);
  cout << *tester->begin() << endl;

  ScopeWaveform* ptrFirstWf = *tester->begin();

  for (int ii=0; ii<5; ii++) {
    // works
    //cout << "element: " << (int)((*tester)[0]->getTrace())[ii] << endl;
    cout << "element: " <<  (int) ptrFirstWf->getTrace()[ii] << endl;
  }
  //cout << "sh->scope()->createHisto()" << endl;

  // test the access to the vector of pointers to scope traces
  // getChannelData() returns a pointer to the vector of pointers.
  /*
  vector<ScopeWaveform*>* tester;
  tester = sh->scope()->getChannelData(1);
  for (vector<ScopeWaveform*>::iterator i=tester->begin();
       i != tester->end(); i++) {
    cout << "hi";
  }
  */

  // deleting vectors used to hold data
  sh->scope()->clearData();

  //sh->scope()->getChannelData(1)

}
