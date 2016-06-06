// Scope.cc
//
#include <iostream>

#include "Scope.hh"

using std::cout;
using std::endl;
//using std::string;
using std::vector;

//const string Scope::ALAZAR_ATS860 = "alazar_ats860";
//const string Scope::INVALID_SCOPE = "invalid";
const TString Scope::ALAZAR_ATS860 = "alazar_ats860";
const TString Scope::INVALID_SCOPE = "invalid";

//const int Scope::SCOPE_SETUP_TEST = 0;
const int Scope::SCOPE_SETUP_TEST(0);

Scope::~Scope() {}

Scope::Scope() {
  cout << GetName() << ": empty constructor" << endl;
  //setNValidTriggers(NULL);
}

Scope::Scope(int boardID, int nchan, int nLevels) {  
  _nChannels = nchan; 
  _nLevels   = nLevels;
  _boardID = boardID;
}

bool
//Scope::isValidType(const string* scopetype) {
Scope::isValidType(const TString* scopetype) {
  if (*scopetype == ALAZAR_ATS860)
    return true;
  else {
    cout << "Invalid scope type -- " << *scopetype << endl;
    return false;
  }
}

char* 
Scope::GetName() { return "Scope"; }

// THESE SHOULD BE IN ScopeData NOT IN Scope
//int 
//Scope::getNValidTriggers() { return _nValidTriggers; }
//void 
//Scope::setNValidTriggers(int nval) { _nValidTriggers=nval; }
//

void 
Scope::setScopeConfig(ScopeConfig* sc) { _scopeConfig=sc; }
ScopeConfig*
Scope::getScopeConfig() { return _scopeConfig; }

int Scope::getNChannels() { return _nChannels; }
// get the voltage resolution in ADC steps 
int Scope::getNLevels() { return _nLevels; }

//
// UNUSED
//

