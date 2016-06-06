// Scope.hh
//
#ifndef SCOPE_HH
#define SCOPE_HH

//class string;
//class vector;
//#include <string>
#include <vector>

#include "TString.h"

#include "AlazarApi.h"
#include "AlazarCmd.h"

class ScopeData;
#include "ScopeConfig.hh"
#include "ScopeTypes.hh"

//#include "ScopeHandler.hh"

//#include "ScopeWaveform.hh"
//#include "ScopeWaveformHeader.hh"
//#include "ScopeRecord.hh"
//#include "TScopeRecord.hh"
//#include "TROOT.h"
//#include "TH1F.h"

//
// base class for computer scopes
// 
class Scope { 

public:

  // Destructor
  virtual ~Scope();
  // Constructors
  Scope();
  // Specify the number of channels and number of levels of voltage resolution
  Scope(int boardID, int nchan, int nlevels);

  // Specific board identifiers
  static const TString ALAZAR_ATS860;
  static const TString INVALID_SCOPE;
  //static const std::string ALAZAR_ATS860;
  //static const std::string INVALID_SCOPE;
  
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

  // To define a pure virtual function (which ensures that the 
  // subclasses will have this function), do
  //virtual type fnName() = 0;
  // verify that a scope is present, get handle of that scope
  // and verify the scope type is as assumed
  virtual int openScope() = 0; 
  virtual ScopeConfig* configureDAQState(DAQ_State daqstate) = 0;
  //virtual int acquireTriggers(ScopeHandler* scope, float duration_ms) = 0;
  //virtual int acquireTriggers(float duration_ms) = 0;
  virtual int acquireTriggers(ScopeData* data, float duration_ms) = 0;
  virtual int readData(ScopeData* data) = 0;
  virtual float getMaxVoltage(U32 voltRangeSpecifier) = 0;
  virtual void setTriggerLevel1Volts(int chanNumber, float level_volts) = 0;
  virtual void setTriggerLevel1(U32 level) = 0;
  virtual U32  getTriggerLevel1() = 0;

  virtual void setTriggerLevel2Volts(int chanNumber, float level_volts) = 0;
  virtual void setTriggerLevel2(U32 level) = 0;
  virtual U32  getTriggerLevel2() = 0;


  ////////////////////////////////////////////////////////
  // trigger set/get
  virtual void setTrigLevel(int itrg, U32 level) = 0;
  virtual U32  getTrigLevel(int itrg) = 0;

  virtual void setTrigEngine(int itrg, ScopeTriggerEngine engine) = 0;
  virtual ScopeTriggerEngine getTrigEngine(int itrg) = 0;

  virtual void setTrigSource(int itrg, ScopeTriggerSource source) = 0;
  virtual ScopeTriggerSource getTrigSource(int itrg) = 0;

  virtual void setTrigSlope(int itrg, ScopeTriggerSlope slope) = 0;
  virtual ScopeTriggerSlope  getTrigSlope(int itrg) = 0;

  virtual void setTrigEngineOperation(ScopeTriggerEngOp operation) = 0;
  virtual ScopeTriggerEngOp  getTrigEngineOperation() = 0;

  virtual void setVoltageRange(int chanNumber, U32 range) = 0;
  virtual U32  getVoltageRange(int chanNumber) = 0;
  virtual void setInputImpedance(int chanNumber, U32 impedance) = 0;
  virtual void setInputCoupling(int chanNumber, ScopeInputCoupling coupling) = 0;
  virtual U32  setTriggerDelay(U32 delay) = 0;
  // end trigger set/get
  ////////////////////////////////////////////////////////

  //virtual float getDtOfClockRate(AU32 sampleRate) = 0;

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

  // THESE SHOULD BE IN ScopeData, NOT IN Scope.
  //int      getNValidTriggers();
  //void     setNValidTriggers(int nval);

  void         setScopeConfig(ScopeConfig* sc);
  ScopeConfig* getScopeConfig();
  int          getNChannels();
  int          getNLevels();  // get the number of ADC steps

  virtual void printScopeSetup() = 0;

  // 
  // UNUSED
  //
  //virtual int initScope() = 0; // pure virtual
  //virtual int setupTrigger(BoardDef bd) = 0;
  //virtual int setupCaptureClock(BoardDef bd) = 0;
  //virtual int setupRecordSize(BoardDef bd) = 0;
  //virtual void makeScopeConfig(const int sc) = 0;
  //virtual void clearData() = 0;
  //virtual std::vector<ScopeWaveform*>* getChannelData(int chNumber) = 0;
  //virtual std::vector<ScopeRecord*>* getRecords(int chNumber) = 0;
  //virtual TScopeRecord* recordOfWaveform(int ii) = 0;
  //virtual TH1F* histoOfWaveform(int iwf) = 0;


protected:
  RETURN_CODE status;

  DAQ_State _state;

  int _boardID;

private:
  //int _nValidTriggers;
  char* GetName();
  ScopeConfig* _scopeConfig;


  // Depends on the hardware.  How can you force this to be set 
  // in the constructor of sub-classes of Scope? 
  // e.g. _nChannels = 2 for Alazar ATS860.
  int _nChannels;  
  int _nLevels;  // number of ADC levels (e.g. 256 for an 8-bit digitizer)

  //
  // UNUSED
  //
  //extern HANDLE h; // board handle  --  may need to extend this for non-alazar boards
  //extern int h;
  //BoardDef _bd;

};


#endif
