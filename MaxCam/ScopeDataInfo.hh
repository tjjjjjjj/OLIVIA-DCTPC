//ScopeDataInfo.hh
//
#ifndef SCOPE_DATA_INFO_HH
#define SCOPE_DATA_INFO_HH

#include "TObject.h"

#include "ScopeTypes.hh"

//
// class to hold scope information about the scope as it
// was configured during data acquisition
//
class ScopeDataInfo : public TObject {

public:

  // Destructor
  ~ScopeDataInfo();

  // Constructors
  ScopeDataInfo();

  // Copy constructor
  //ScopeDataInfo(const ScopeDataInfo &other);

  // set
  void setSamplingRate(float sr);
  void setVoltageRange(float vmin, float vmax);
  void setVoltageMin(float vmin);
  void setVoltageMax(float vmax);
  void setInputCoupling(ScopeInputCoupling ic);
  void setInputImpedance(float imped);
  void setChannelId(int chanId);
  void setNTriggers(int nt);
  void setScopeNumber(int isco); 

  void setVoltageStep(float dv);   // voltage resolution
  void setTriggerLevel(float tl);  // trigger level in *volts*
  void setTriggerSlope(ScopeTriggerSlope sts);  // trigger slope
  //void setTriggerCause(unsigned int cause); // record trigger cause

  // Future
  // setScopeUniqueId();  // needs new SDK (5.x)
  // setClockType(ScopeClockSource ct); ??? internal/external?
  // setClockEdgeType(); CLOCK_EDGE_RISING or CLOCK_EDGE_FALLING
  // setTriggerType()???
  // 

  // get
  float getSamplingRate() const;
  float getVoltageMin() const;
  float getVoltageMax() const;
  //float getChannelId();
  int   getChannelId() const;

  float getInputImpedance() const;
  int   getScopeNumber() const; 
  int   getNTriggers() const;
  ScopeInputCoupling getInputCoupling() const;
  float getVoltageStep() const;
  float getTriggerLevel() const;
  ScopeTriggerSlope getTriggerSlope() const;
  //unsigned int getTriggerCause();   

  const char* GetName() const;

private:
  float _samplingRate;
  float _voltageMin;
  float _voltageMax;
  float _voltageStep;
  int   _nTriggers;
  float _triggerLevelVolts;
  ScopeTriggerSlope _triggerSlope;
  ScopeInputCoupling _inputCoupling;
  float _inputImpedance;
  int _channelId;
  int _scopeNumber;
  //unsigned int _triggerCause; //

ClassDef(ScopeDataInfo, 2)
};

#endif
