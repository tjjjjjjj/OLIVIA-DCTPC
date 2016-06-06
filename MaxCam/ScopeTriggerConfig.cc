//ScopeTriggerConfig.cc
//
#include <iostream>

#include "ScopeTriggerConfig.hh"

ClassImp(ScopeTriggerConfig)

ScopeTriggerConfig::~ScopeTriggerConfig() {}
ScopeTriggerConfig::ScopeTriggerConfig() {}
ScopeTriggerConfig::ScopeTriggerConfig(AU32 trgOp, 
				       AU32 trgEn1, AU32 trgSrc1, 
				       AU32 trgSlope1, AU32 trgLevel1,
				       AU32 trgEn2, AU32 trgSrc2, 
				       AU32 trgSlope2, AU32 trgLevel2) {
  _trgOp     = trgOp;
  _trgEn1    = trgEn1;
  _trgSrc1   = trgSrc1;
  _trgSlope1 = trgSlope1;
  _trgLevel1 = trgLevel1;
  _trgEn2    = trgEn2;
  _trgSrc2   = trgSrc2;
  _trgSlope2 = trgSlope2;
  _trgLevel2 = trgLevel2;
}

AU32 ScopeTriggerConfig::trigOp() { return _trgOp; }
AU32 ScopeTriggerConfig::trigEn1() { return _trgEn1; }
AU32 ScopeTriggerConfig::trigEn2() { return _trgEn2; }
AU32 ScopeTriggerConfig::trigSrc1() { return _trgSrc1; }
AU32 ScopeTriggerConfig::trigSrc2() { return _trgSrc2; }
AU32 ScopeTriggerConfig::trigSlope1() { return _trgSlope1; }
AU32 ScopeTriggerConfig::trigSlope2() { return _trgSlope2; }
AU32 ScopeTriggerConfig::trigLevel1() { return _trgLevel1; }
AU32 ScopeTriggerConfig::trigLevel2() { return _trgLevel2; }

// setters
void ScopeTriggerConfig::setTrigOp(AU32 to) { _trgOp = to; }

void ScopeTriggerConfig::setTrigLevel1(AU32 tl) { _trgLevel1 = tl; }
void ScopeTriggerConfig::setTrigEngine1(AU32 te) { _trgEn1 = te; }
void ScopeTriggerConfig::setTrigSource1(AU32 ts) { _trgSrc1 = ts; }
void ScopeTriggerConfig::setTrigSlope1(AU32 ts) { _trgSlope1 = ts; }

void ScopeTriggerConfig::setTrigLevel2(AU32 tl) { _trgLevel2 = tl; }
void ScopeTriggerConfig::setTrigEngine2(AU32 te) { _trgEn2 = te; }
void ScopeTriggerConfig::setTrigSource2(AU32 ts) { _trgSrc2 = ts; }
void ScopeTriggerConfig::setTrigSlope2(AU32 ts) { _trgSlope2 = ts; }

