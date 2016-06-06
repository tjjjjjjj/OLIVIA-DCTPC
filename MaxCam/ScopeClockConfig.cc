//ScopeClockConfig.cc
//
#include <iostream>

#include "ScopeClockConfig.hh"

ClassImp(ScopeClockConfig)

using std::cout;
using std::endl;
using std::string;

ScopeClockConfig::~ScopeClockConfig() {}
ScopeClockConfig::ScopeClockConfig() {}
ScopeClockConfig::ScopeClockConfig(AU32 cs, AU32 cr, AU32 ce, AU32 cd) {
    setClockSource(cs);
    setClockRate(cr);
    setClockEdge(ce);
    setClockDecimation(cd);
  }

void ScopeClockConfig::setClockSource(AU32 cs) { _clockSource = cs; }
void ScopeClockConfig::setClockRate(AU32 cr) { _clockRate = cr; }
void ScopeClockConfig::setClockEdge(AU32 ce) { _clockEdge = ce; }
void ScopeClockConfig::setClockDecimation(AU32 cd) { _clockDecimation = cd; }

AU32 ScopeClockConfig::source() { return _clockSource; } // is this used?
AU32 ScopeClockConfig::getClockSource() { return _clockSource; }
AU32 ScopeClockConfig::getClockRate()   { return _clockRate; }
AU32 ScopeClockConfig::getClockEdge()   { return _clockEdge; }
AU32 ScopeClockConfig::getClockDecimation() { return _clockDecimation; }
