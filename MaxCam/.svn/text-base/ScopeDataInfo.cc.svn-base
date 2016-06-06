//ScopeDataInfo.cc
//
#include <iostream>
#include "ScopeDataInfo.hh"

ClassImp(ScopeDataInfo)

using std::cout;
using std::endl;

ScopeDataInfo::~ScopeDataInfo() {
  //cout << GetName() << ": destructor" << endl;
}

ScopeDataInfo::ScopeDataInfo() {}

//// copy constructor
//ScopeDataInfo::ScopeDataInfo(const ScopeDataInfo &other) {
//  setSamplingRate(other._samplingRate);
//  //_samplingRate(other._samplingRate);
//}

// set
void ScopeDataInfo::setSamplingRate(float sr) { _samplingRate = sr; };
void ScopeDataInfo::setVoltageMin(float vmin) { _voltageMin = vmin; }
void ScopeDataInfo::setVoltageMax(float vmax) { _voltageMax = vmax; }
void ScopeDataInfo::setChannelId(int chanId) { _channelId = chanId; }
void ScopeDataInfo::setInputCoupling(ScopeInputCoupling ic) { _inputCoupling = ic; }
void ScopeDataInfo::setInputImpedance(float imped) { _inputImpedance = imped; }
void ScopeDataInfo::setNTriggers(int nt) { _nTriggers = nt; }
void ScopeDataInfo::setScopeNumber(int isco) { _scopeNumber = isco; }
void ScopeDataInfo::setVoltageRange(float vmin, float vmax) {
  setVoltageMin(vmin);
  setVoltageMax(vmax);
}
void ScopeDataInfo::setVoltageStep(float vstep) { _voltageStep = vstep; }
void ScopeDataInfo::setTriggerLevel(float tl) { _triggerLevelVolts = tl; }
void ScopeDataInfo::setTriggerSlope(ScopeTriggerSlope sts) { _triggerSlope = sts; }

// get
int   ScopeDataInfo::getChannelId() const{ return _channelId; }
float ScopeDataInfo::getSamplingRate() const{ return _samplingRate; }
float ScopeDataInfo::getVoltageMin() const{ return _voltageMin; }
float ScopeDataInfo::getVoltageMax() const{ return _voltageMax; }
float ScopeDataInfo::getInputImpedance() const{ return _inputImpedance; }
int   ScopeDataInfo::getScopeNumber() const{ return _scopeNumber; }
int   ScopeDataInfo::getNTriggers() const{ return _nTriggers; }
ScopeInputCoupling ScopeDataInfo::getInputCoupling() const{ return _inputCoupling; }
float ScopeDataInfo::getVoltageStep() const{ return _voltageStep; }
float ScopeDataInfo::getTriggerLevel() const{ return _triggerLevelVolts; }
ScopeTriggerSlope ScopeDataInfo::getTriggerSlope() const{ return _triggerSlope; }

// 
// PRIVATE
//
const char* ScopeDataInfo::GetName() const { return "ScopeDataInfo"; }

