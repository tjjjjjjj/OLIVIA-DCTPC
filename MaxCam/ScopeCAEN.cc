// ScopeCAEN.cc
//
#include <iostream>

#include "ScopeCAEN.hh"

using std::cout;
using std::endl;
using std::vector;

const TString ScopeCAEN::CAENv1720 = "caenv1720";
const TString ScopeCAEN::INVALID_SCOPE = "invalid";

const int ScopeCAEN::SCOPE_SETUP_TEST(0);

ScopeCAEN::~ScopeCAEN() {}

ScopeCAEN::ScopeCAEN() {
  cout << GetName() << ": empty constructor" << endl;
}

ScopeCAEN::ScopeCAEN(int boardID, int nchan, int nLevels) {  
  _nChannels = nchan; 
  _nLevels   = nLevels;
  _boardID = boardID;
}

bool
ScopeCAEN::isValidType(const TString* scopetype) {
  if (*scopetype == CAENv1720)
    return true;
  else {
    cout << "Invalid scope type -- " << *scopetype << endl;
    return false;
  }
}

char* 
ScopeCAEN::GetName() { return "CAENv1720"; }

void 
ScopeCAEN::setScopeConfig(ScopeConfig* sc) { _scopeConfig=sc; }
ScopeConfig*
ScopeCAEN::getScopeConfig() { return _scopeConfig; }

int ScopeCAEN::getNChannels() { return _nChannels; }
// get the voltage resolution in ADC steps 
int ScopeCAEN::getNLevels() { return _nLevels; }

//
// UNUSED
//

