// ScopeAlazarATS860.hh
//
#ifndef SCOPE_ALAZAR_ATS860_HH
#define SCOPE_ALAZAR_ATS860_HH

//#include <string>
//#include <vector>

#include "TROOT.h"
#include "TString.h"

#include "AlazarApi.h"
#include "AlazarCmd.h"

#include "Scope.hh"
//#include "ScopeTypes.hh"

typedef long long INT64;

//
// Alazar ATS860 Class (PCI Digitizer card -- an oscilloscope for your computer)
// 
// 
class ScopeAlazarATS860 : public Scope { 

public:

  // Destructor
  ~ScopeAlazarATS860();

  // Constructor
  ScopeAlazarATS860(int boardID);

  int openScope();
  int acquireTriggers(ScopeData* data, float duration_ms);
  int readData(ScopeData* data);
  float getMaxVoltage(U32 voltRangeSpecifier);

  // returns time span of each sample in seconds
  //float getDtOfClockRate(AU32 sampleRate);  

  //////////////
  // used by ScopeDataInfo class
  //
  float getSamplingRate();
  float getVoltageMin(int chanNum);
  float getVoltageMax(int chanNum);
  float getInputImpedance(int chanNum);
  ScopeInputCoupling getInputCoupling(int chanNum);
  float getVoltageStep(int chanNum);
  float getTriggerLevel(int chanNum);
  ScopeTriggerSlope getTriggerSlope(int chanNum);
  unsigned int setTriggerSlope(unsigned int slope);
  //

  ScopeConfig* configureDAQState(Scope::DAQ_State daqstate);

private:
  const static int ACQ_MODE_TRADITIONAL;
  const static int ACQ_MODE_AUTODMA;
  const static int N_RECORDS_PER_BUFFER;
  static U8 const ALAZAR_ATS860_CHANNELS[];

  int applyScopeConfig(ScopeConfig *config);


  HANDLE _h;
  U32 UseHeader;
  U32 channelMode;
  U32 RecsPerBuffer;
  BoardDef     getBoardDef(Scope::DAQ_State daqstate);
  //ScopeConfig* makeScopeConfig(BoardDef* bd, std::string* chALabel, std::string* chBLabel);
  ScopeConfig* makeScopeConfig(BoardDef* bd, TString* chALabel, TString* chBLabel);

  TString _chALabel, _chBLabel;

  // buffers
  U8* _buffers[2];
  void allocateMemoryForBuffers(ScopeConfig *sc);

  //
  int configureHardwareForDAQ(ScopeConfig* sc);
  int setRecordSize(ScopeConfig* sc);
  int setRecordCount(ScopeConfig* sc);
  int setCaptureClock(ScopeConfig* sc);
  int setChannel(U8 chan, ScopeConfig* sc);
  int setTrigger(ScopeConfig* sc);
  int computeBytesPerRecord(BoardDef* bd);

  U32   levelOfVolts(int chanNumber, float level_volts);
  float voltsOfLevel(int chanNumber, U32   level);
  
  void setTriggerLevel1Volts(int chanNumber, float level_volts);
  void setTriggerLevel2Volts(int chanNumber, float level_volts);
  void setTriggerLevel1(U32 level);
  void setTriggerLevel2(U32 level);
  U32  getTriggerLevel1();
  U32  getTriggerLevel2();
  void setVoltageRange(int chanNumber, U32 range);
  U32  getVoltageRange(int chanNumber);
  void setInputImpedance(int chanNumber, U32 impedance);
  void setInputCoupling(int chanNumber, ScopeInputCoupling coupling);
  U32  setTriggerDelay(U32 delay);


  ///////////////////////////////////////////////////////
  // Trigger set/get methods
  // These supercede the older ones...
  // jbattat 8/24/2010
  void setTrigLevel(int itrg, U32 level);
  U32  getTrigLevel(int itrg);

  void setTrigEngine(int itrg, ScopeTriggerEngine engine);
  ScopeTriggerEngine getTrigEngine(int itrg);
  
  void setTrigSource(int itrg, ScopeTriggerSource source);
  ScopeTriggerSource  getTrigSource(int itrg);
  
  void setTrigSlope(int itrg, ScopeTriggerSlope slope);
  ScopeTriggerSlope  getTrigSlope(int itrg);
  
  void setTrigEngineOperation(ScopeTriggerEngOp operation);
  ScopeTriggerEngOp getTrigEngineOperation();
  // end of set/get trig methods
  ///////////////////////////////////////////////////////

  void   verifyBoardsExist();
  void   verifyBoardType(HANDLE hh);
  void   printDriverVersion();
  void   printSDKVersion();
  void   printScopeSerialNumber(HANDLE hh);
  void   printScopeMemorySize(HANDLE hh);
  void   printScopeChannelInfo(HANDLE hh);
  void   printScopeSetup();
  int    apiError(RETURN_CODE status);
  void   setScopeHandle(HANDLE h);
  void   freeMemory();
  HANDLE getScopeHandle();

  int acquireTriggersTraditional(float duration_ms);
  //int acquireTriggersAutoDMA();

  double timestampToSeconds(U32 uTimeStampHighPart, U32 uTimeStampLowPart, double dSamplesPerSecond);

  void setDataAcqMode(int mode);
  int _acqMode;

  char* GetName();

  // functions to convert between Alazar specific values and
  // "generic" Scope Type values
  // See, e.g. ScopeTypes.hh
  ScopeTriggerSource AlazarToGenericSource(AU32 source);
  AU32 GenericToAlazarSource(ScopeTriggerSource source);
  ScopeTriggerSlope AlazarToGenericSlope(AU32 slope);
  AU32 GenericToAlazarSlope(ScopeTriggerSlope slope);
  ScopeTriggerEngine AlazarToGenericEngine(AU32 engine);
  AU32 GenericToAlazarEngine(ScopeTriggerEngine engine);
  ScopeTriggerEngOp AlazarToGenericEngOp(AU32 engOp);
  AU32 GenericToAlazarEngOp(ScopeTriggerEngOp engOp);

  //
  // below here is unused i think
  //
  int processBufferData(U8* Data);

 

};

#endif
