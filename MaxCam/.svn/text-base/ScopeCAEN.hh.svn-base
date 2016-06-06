// Scope.hh
//
#ifndef SCOPE_CAEN_HH
#define SCOPE_CAEN_HH

//class string;
//class vector;
//#include <string>
#include <vector>

#include "TString.h"

#ifdef SCOPE_CAEN
#include <CAENComm.h>
#include <CAENDigitizer.h>
#else
#include "AlazarApi.h"
#include "AlazarCmd.h"
#endif


class ScopeData;
#include "ScopeConfig.hh"
#include "ScopeTypes.hh"

//
// base class for computer scopes
// 
class ScopeCAEN { 

public:

  // Destructor
  virtual ~ScopeCAEN();
  // Constructors
  ScopeCAEN();
  // Specify the number of channels and number of levels of voltage resolution
  ScopeCAEN(int boardID, int nchan, int nlevels);

  // Specific board identifiers
  static const TString CAENv1720;
  static const TString INVALID_SCOPE;
   
  enum DAQ_State { DAQ_NONE, // no scope 
		   DAQ_CHARGE, // legacy flag
		   DAQ_TEST, // legacy flag
		   DAQ_CV1, // testing for z-coordinate
		   DAQ_TGA_ROA, // trigger A, read A 
		   DAQ_TGA_ROAB=DAQ_CHARGE, // trigger A, read A&B
		   DAQ_TGAB_ROAB, // trigger A read A, trigger B read B
		   DAQ_TGEXT_ROAB // external trigger, read A&B
  };
  void setDAQState(DAQ_State state) { _state=state; }
  DAQ_State getDAQState() { return _state; }

  const static int SCOPE_SETUP_TEST;

  enum ChanID {CHAN_A, CHAN_B, CHAN_C, CHAN_D};

  // and verify the scope type is as assumed
  virtual int openScope() = 0; 
  virtual ScopeConfig* configureDAQState(DAQ_State daqstate) = 0;
  virtual int acquireTriggers(ScopeData* data, float duration_ms) = 0;
  virtual int readData(ScopeData* data) = 0;
  virtual void setTriggerLevel1Volts(int chanNumber, float level_volts) = 0;

  virtual void setTrigEngine(int itrg, ScopeTriggerEngine engine) = 0;
  virtual ScopeTriggerEngine getTrigEngine(int itrg) = 0;

  virtual void setTrigSource(int itrg, ScopeTriggerSource source) = 0;
  virtual ScopeTriggerSource getTrigSource(int itrg) = 0;

  virtual void setTrigSlope(int itrg, ScopeTriggerSlope slope) = 0;
  virtual ScopeTriggerSlope  getTrigSlope(int itrg) = 0;

  virtual void setTrigEngineOperation(ScopeTriggerEngOp operation) = 0;
  virtual ScopeTriggerEngOp  getTrigEngineOperation() = 0;

  //static bool isValidType(const std::string* scopetype);
  static bool isValidType(const TString* scopetype);

  //////////////
  // called by ScopeDataInfo class
  virtual float getSamplingRate() = 0;
  virtual float getVoltageMin(int chanNum) = 0;
  virtual float getVoltageMax(int chanNum) = 0;
  virtual float getInputImpedance(int chanNum) = 0;
  virtual ScopeInputCoupling getInputCoupling(int chanNum) = 0;
  virtual float getVoltageStep(int chanNum) = 0;
  virtual float getTriggerLevel(int chanNum) = 0;
  virtual ScopeTriggerSlope getTriggerSlope(int chanNum) = 0;
  virtual unsigned int setTriggerSlope(unsigned int slope) = 0;
  /////////////

  void         setScopeConfig(ScopeConfig* sc);
  ScopeConfig* getScopeConfig();
  int          getNChannels();
  int          getNLevels();  // get the number of ADC steps

  virtual void printScopeSetup() = 0;

protected:
  //RETURN_CODE status;

  DAQ_State _state;

  int _boardID;

private:
  char* GetName();
  ScopeConfig* _scopeConfig;


  // Depends on the hardware.  How can you force this to be set 
  // in the constructor of sub-classes of Scope? 
  // e.g. _nChannels = 2 for Alazar ATS860.
  int _nChannels;  
  int _nLevels;  // number of ADC levels (e.g. 256 for an 8-bit digitizer)

};


#endif
