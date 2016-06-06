// ScopeHandler.hh
//
#ifndef SCOPE_HANDLER_HH
#define SCOPE_HANDLER_HH

//#include <string>
#include <vector>

#include "TString.h"
#include <mysql/mysql.h>
#include "AlazarApi.h"
#include "AlazarCmd.h"

#include "ScopeTypes.hh"
//Alazar types:
// U8  is unsigned char
// U32 is unsigned int


#include "Scope.hh"
//#include "ScopeData.hh"
class ScopeData;
class ScopeDataInfo;
class ScopeWaveformData;
class ScopeConfig;

//class ScopeAlazarATS860;
//#include "ScopeAlazarATS860.hh"

//
// Class to run an oscilloscope
//
class ScopeHandler {

public:

  // Destructor
  ~ScopeHandler();
  // Constructor
  ScopeHandler();
  ScopeHandler(const TString* scopeType);
  //ScopeHandler(const std::string* scopeType);

  void addScope(const TString* scopeType);
  //void addScope(const std::string* scopeType);
  void openScope(int ii);
  void openScope();

  ScopeConfig* configureDAQState(int ii=0, Scope::DAQ_State daqstate=Scope::DAQ_CHARGE);

  ///////////////////////////////////////////////////////////
  // Trigger set/get using Alazar-specific types....
  void setTrigLevel(int isco, int itrg, unsigned int level);
  void setTrigSlope(int isco, int itrg, ScopeTriggerSlope slope);
  void setTrigSource(int isco, int itrg, ScopeTriggerSource source);
  void setTrigEngine(int isco, int itrg, ScopeTriggerEngine engine);
  void setTrigEngineOperation(int isco, ScopeTriggerEngOp operation);

  unsigned int getTrigLevel(int isco, int itrg);
  ScopeTriggerSlope getTrigSlope(int isco, int itrg);
  ScopeTriggerSource getTrigSource(int isco, int itrg);
  ScopeTriggerEngine getTrigEngine(int isco, int itrg);
  ScopeTriggerEngOp getTrigEngineOperation(int isco);

  void setTriggerLevel1Volts(int ii, int chanNumber, float level_volts);
  unsigned int setTriggerDelay( int iScope, unsigned int delay);
  unsigned int setTriggerSlope( int iSlope, unsigned int slope);

  // deprecated... see setTriggerLevel(int itrg, unsigned int level)
  void setTriggerLevel1(int ii, unsigned int level);
  unsigned int  getTriggerLevel1(int ii);
  void setTriggerLevel2(int ii, unsigned int level);
  unsigned int  getTriggerLevel2(int ii);
  // end trigger set/get

  void setVoltageRange(int scopeNumber, int chanNumber, unsigned int range);
  unsigned int  getVoltageRange(int scopeNumber, int chanNumber);

  void setInputImpedance(int scopeNumber, int chanNumber, unsigned int impedance);
  void setInputCoupling(int scopeNumber, int chanNumber, ScopeInputCoupling coupling);

  //void setTriggerCause(unsigned int cause);
  //unsigned int getTriggerCause();

  const TString* getScopeType(int ii=0);
  //const std::string* getScopeType(int ii=0);

  ScopeConfig*  getScopeConfig(int ii=0);

  /////////////////
  // called by ScopeDataInfo class
  float getSamplingRate(int isco);
  float getVoltageMin(int isco, int chanNum);
  float getVoltageMax(int isco, int chanNum);
  float getInputImpedance(int isco, int chanNum);
  ScopeInputCoupling getInputCoupling(int isco, int chanNum);
  float getVoltageStep(int isco, int chanNum);
  float getTriggerLevel(int isco, int chanNum); 
  ScopeTriggerSlope getTriggerSlope(int isco, int chanNum);
  //

  int getNChannels(int isco);
  int getNLevels(int isco);

  int  acquireTriggers(int iscope, float duration_ms=1000.0);
  int  readData(int ii);
  void addScopeData(int ii, ScopeConfig* sc);
  void addScopeDataInfo(int ii);

  ScopeData* data(int ii=0);
  ScopeDataInfo* dataInfo(int ii=0);
  
  int getNValidTriggers(int ii=0);
  std::vector<ScopeWaveformData*> getWaveforms(int ii=0);

  float getMaxVoltage(int ii, unsigned int voltRangeSpecifier);
  //float getDtOfClockRate(int ii, AU32 sampleRate);

  void clearData(int ii=0);

  unsigned int nScope() { return _scopeList.size(); }

  enum ScopeReadoutScheme {
    TRIGGER_CHA_READOUT_CHA, // trigger on A, read A only
    TRIGGER_CHA_READOUT_CHA_CHB, // trigger on A, read A&B
    TRIGGER_CHA_CHB_READOUT_CHA_CHB, // trigger A and B separately
    TRIGGER_EXT_READOUT_CHA, // external trigger, read A
    TRIGGER_EXT_READOUT_CHA_CHB // external trigger, read A&B
  };

  void saveScopeInfoToDB(const char *fname, MYSQL * db_handle = NULL);

private:
  // return pointer to the scope 
  // currently just returns the first scope in the list
  Scope* scope(int ii=0);

  void validateScopeType(const TString* scopetype);
  void setScopeType(const TString* scopeType);
  //void validateScopeType(const std::string* scopetype);
  //void setScopeType(const std::string* scopeType);
  
  // vector of pointers to instances of the Scope class
  // or subclasses of Scope (e.g. ScopeAlazarATS860)
  //std::vector<const std::string*> _scopeTypes;
  std::vector<Scope*>         _scopeList;
  std::vector<const TString*> _scopeTypes;
  std::vector<ScopeData*>     _scopeDatasets;
  std::vector<ScopeDataInfo*> _scopeDataInfo;

  char* GetName();
};


#endif
