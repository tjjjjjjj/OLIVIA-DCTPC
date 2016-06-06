//ScopeChannelConfig.cc
//
#include "ScopeChannelConfig.hh"

#include <iostream>

using std::cout;
using std::endl;
//using std::string;

ClassImp(ScopeChannelConfig)

ScopeChannelConfig::~ScopeChannelConfig() {
  cout << GetName() << ": destructor" << endl;
  // do not "delete _label" because "_label" in multiple
  // instances of this class could be pointing to the same
  // address location...
  //delete _label;  
}

ScopeChannelConfig::ScopeChannelConfig() {}
ScopeChannelConfig::ScopeChannelConfig(AU8  chanNum, AU32 coupling, 
				       AU32 vRange,  AU32 impedance, 
				       TString* label) {
				       //string* label) {
  _channelNumber = chanNum;
  _coupling      = coupling;
  _voltageRange  = vRange;
  _impedance     = impedance;
  _label         = label;  // e.g. "charge"
}

//string* ScopeChannelConfig::label() { return _label; }
TString* ScopeChannelConfig::label() { return _label; }

AU8  ScopeChannelConfig::getChannelNumber()  { return _channelNumber; }
AU32 ScopeChannelConfig::getCoupling()       { return _coupling; }
AU32 ScopeChannelConfig::getVoltageRange()   { return _voltageRange; }
AU32 ScopeChannelConfig::getInputImpedance() { return _impedance; }

void ScopeChannelConfig::setVoltageRange(AU32 range) { _voltageRange = range; }
void ScopeChannelConfig::setInputCoupling(AU32 coupling) {_coupling = coupling; }
void ScopeChannelConfig::setInputImpedance(AU32 impedance) {_impedance = impedance; }

// 
// PRIVATE
//
char* ScopeChannelConfig::GetName() { return "ScopeChannelConfig"; }
