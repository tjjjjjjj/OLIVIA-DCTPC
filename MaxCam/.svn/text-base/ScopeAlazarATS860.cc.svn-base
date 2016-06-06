

// ScopeAlazarATS860.cc
//
#include <iostream>
#include <assert.h>
#include <stdlib.h>
//#include <ctime>
#include "ScopeAlazarATS860.hh"
#include "ScopeHandler.hh"
#include "ScopeWaveform.hh"
#include "ScopeData.hh"
//#include "TString.h"
#include "TStopwatch.h"

using std::cout;
using std::endl;
//using std::string;
using std::ofstream;
using std::ios;
using std::flush;
using std::vector;

const int ScopeAlazarATS860::ACQ_MODE_TRADITIONAL = 0;
const int ScopeAlazarATS860::ACQ_MODE_AUTODMA = 1;
const int ScopeAlazarATS860::N_RECORDS_PER_BUFFER = 25;
//const int ScopeAlazarATS860::N_RECORDS_PER_BUFFER = 2;

U8 const ScopeAlazarATS860::ALAZAR_ATS860_CHANNELS[2] = {CHANNEL_A, CHANNEL_B};

//
// PUBLIC Methods
//

ScopeAlazarATS860::~ScopeAlazarATS860() {
  free(_buffers[0]);
  free(_buffers[1]);
};

ScopeAlazarATS860::ScopeAlazarATS860(int boardID) : Scope(boardID, 2, 256) {
  cout << GetName() << ": constructor" << endl;

  // Catch API errors
  status = ApiSuccess;

  setDataAcqMode(ACQ_MODE_TRADITIONAL);
  //setDataAcqMode(ACQ_MODE_AUTODMA);

}

void
ScopeAlazarATS860::allocateMemoryForBuffers(ScopeConfig* sc) {

  if (_acqMode == ACQ_MODE_TRADITIONAL) {
    int bytesPerRecord = sc->board()->getBytesPerRecord();
    int nChannels = 2;
    for (int ii=0; ii<nChannels; ii++) {
      _buffers[ii] = (U8 *) malloc(bytesPerRecord);
    }
  } else if (_acqMode == ACQ_MODE_AUTODMA) {
    cout << "ACQ_MODE_AUTODMA not yet implemented" << endl;
    cout << "No buffer allocation done yet." << endl;
    cout << "see " << GetName() << ".allocateMemoryForBuffers()" << endl;
    assert(!"allocateMemoryForBuffers error:  ACQ_MODE_AUTODMA not implemented");
  } else {
    assert(!"allocateMemoryForBuffers error:  unrecognized _acqMode");
  }

}

void
ScopeAlazarATS860::setDataAcqMode(int mode) {
  // AUTODMA mode can collect data and write data
  // with no dead time.  Faster than TRADITIONAL mode (serial)
  // but not supported in Linux yet...

  _acqMode = mode;

  if (_acqMode == ACQ_MODE_TRADITIONAL) {
    cout << "Data acquisition mode = ACQ_MODE_TRADITIONAL = " 
	 << ACQ_MODE_TRADITIONAL << endl;
  } else if (_acqMode == ACQ_MODE_AUTODMA) {
    // Set variables that are only used in AUTODMA mode
    cout << "Data acquisition mode = ACQ_MODE_AUTODMA = " 
	 << ACQ_MODE_AUTODMA << endl;

    // get header information as well (!= 0 means yes, use headers)
    UseHeader = 1;

    channelMode = CHANNEL_A | CHANNEL_B;  // dual channel mode
    // channelMode = CHANNEL_A  // channel A only
    // channelMode = CHANNEL_B  // channel B only
  }

}

int
ScopeAlazarATS860::openScope() {
  // establish communication with scope
  // and save scope handle for later use

  cout << GetName() << ":openScope()" << endl;

  // Useful for debugging
  printDriverVersion();
  printSDKVersion();

  verifyBoardsExist();

  // get the board handle (for now, assume just one board)
  U32 systemId=1;
  HANDLE hh = AlazarGetBoardBySystemID(systemId, _boardID+1); // assume one system, multiple boards
  cout <<"Trying to get handle using ID = " << (_boardID+1) << "  get = " << hh << endl;

  verifyBoardType(hh);
  printScopeSerialNumber(hh);  // cannot use until SDK 5.X is available
  printScopeMemorySize(hh);
  printScopeChannelInfo(hh);

  // blink the LED (just for fun)
  AlazarSetLED(hh, LED_ON);
  AlazarSetLED(hh, LED_OFF);

  // save the scope handle for later use
  setScopeHandle( hh );

  return 0;
}


int 
ScopeAlazarATS860::applyScopeConfig(ScopeConfig *scopeConfig) {
  // Applies requested configuration.

  // set configuration pointer
  setScopeConfig(scopeConfig);
  // and actually configure the hardware, too
  configureHardwareForDAQ(scopeConfig);
  // and book the buffers for data transfer
  allocateMemoryForBuffers(scopeConfig);

  return 0;
}

ScopeConfig*
ScopeAlazarATS860::configureDAQState(Scope::DAQ_State daqstate) {
  // Make scope configuration for given daq state.

  //cout << GetName() << ": configureDAQState()" << endl;
  
  _state=daqstate;

  ScopeConfig* scopeConfig;

  BoardDef bd;
  if (daqstate == Scope::DAQ_CHARGE) {
    bd = getBoardDef(daqstate);
    _chALabel = "chargeA";
    _chBLabel = "chargeB";
    //scopeConfig = makeScopeConfig(&bd, &_chALabel, &_chBLabel);
    //setScopeConfig(scopeConfig);
    // and actually configure the hardware, too
    //configureHardwareForDAQ(scopeConfig);
  } else if (daqstate == Scope::DAQ_TEST) {
    cout << "******** WARNING::: Scope is set up DAQ_TEST *****" << endl;
    bd = getBoardDef(daqstate);
    _chALabel = "charge";
    _chBLabel = "junk";
  } else if (daqstate == Scope::DAQ_CV1) {
    cout << "Scope is set up DAQ_CV1" << endl;
    bd = getBoardDef(daqstate);
    _chALabel = "charge";
    _chBLabel = "junk";
  } else {
    cout << GetName() << ".configureDAQState() -- ERROR\n" 
	 << "Unknown Scope::DAQ_State --> " << daqstate << endl;
    assert(!"configureDAQState error:  unrecognized DAQ_State");
  }

  scopeConfig = makeScopeConfig(&bd, &_chALabel, &_chBLabel);
  applyScopeConfig(scopeConfig);

  return scopeConfig;
}

int 
ScopeAlazarATS860::acquireTriggers(ScopeData* data, float duration_ms) {
  // duration -- length of time over which to acquire triggers (seconds)

  //cout << GetName() << ".acquireTriggers()" << endl;

  if (_acqMode == ACQ_MODE_TRADITIONAL) {
    return acquireTriggersTraditional(duration_ms);
    //return acquireTriggersTraditional(data, duration_ms);
    //return acquireTriggersTraditional(scope, duration_ms);
  } else if (_acqMode == ACQ_MODE_AUTODMA) {
    //cout << "you asked for AUTODMA acquisition mode\n"
    //     << "This has not yet been implemented yet...\n"
    //	 << "returning -1" << endl;
    assert(!"AUTODMA mode not yet implemented");
    return -1;
  } else {
    cout << "unknown _acqMode value: " << _acqMode
    	 << "returning -1" << endl;
    assert(!"ERROR in acquireTriggers(): unknown _acqMode value");
    return -1;
  }

  return 0;
}



// provided by Nigel Ferdinand of AlazarTech
double 
ScopeAlazarATS860::timestampToSeconds (
		    U32    uTimeStampHighPart, // from AlazarGetTriggerAddress
		    U32    uTimeStampLowPart, // from AlazarGetTriggerAddress
		    double dSamplesPerSecond// sample rate 
		    ) 
{
  INT64 timestamp_counts;
  double timestamp_sec;
  const int multiplier = 4; // constant for the ATS860

  //cout << "uTimeStampHighPart = " << (long)uTimeStampHighPart << endl;
  //cout << "uTimeStampLowPart = " << (long)uTimeStampLowPart << endl;

  timestamp_counts  = ((INT64) uTimeStampHighPart) << 8;
  timestamp_counts |= (INT64) (uTimeStampLowPart & 0x0ff);

  if (dSamplesPerSecond > 0.)
    timestamp_sec = multiplier * timestamp_counts /
      dSamplesPerSecond;
  else
    timestamp_sec = 0.;

  return timestamp_sec;
}

int
ScopeAlazarATS860::acquireTriggersTraditional(float duration_ms) {
  // This call should be done only for master board
  // duration_ms -- how long to acquire triggers for (milliseconds)
  cout << GetName() << ": Acquiring triggers for " << duration_ms
       << " millisecond(s)" << endl;
  
  HANDLE hh = getScopeHandle();

  //
  // The idea here is to set up board to acquire some number of triggers 
  // (>> 1) and then to start taking triggers, pause for the "exposure time"
  // and then stop the acquisition
  // All data will be stored in the on-board memory (PCI).  
  // The Alazar ATS860 board has 16 MB per channel of acquisition memory
  // (can be upgraded to 256 MB/channel)
  // Then, you can query the board to see how many triggers it actually acquired
  // and then you can transfer the data from those triggers
  // 
  // get timestamp right before capturing begins
  //time_t scopeStartTime = time(NULL);
  status = AlazarStartCapture(hh);
  usleep(int(duration_ms*1000));
  status = AlazarAbortCapture(hh);

  return 0;
}

int
ScopeAlazarATS860::readData(ScopeData* data) {

  HANDLE hh = getScopeHandle();

  HANDLE hmaster = AlazarGetBoardBySystemID(1, 1);


  ScopeConfig* sc = getScopeConfig();
  
  //sc->printall();

  int elementSize  = 1; // =1 for 8-bit digitizers, 
  //                    // =2 for 12,a14 and 16 bit.

  int chAId = 0;
  int chBId = 1;
  ScopeDataChannel* dataChanA = data->dataChan(chAId);
  ScopeDataChannel* dataChanB = data->dataChan(chBId);
  int wfAId=0;
  int wfBId=0;

  // Query the board to find out how many 
  // triggers it actually acquired (forced or otherwise)
  U32 uRecords;
  status = AlazarQueryCapability(hh, NUMBER_OF_RECORDS, 0, &uRecords);
  //cout << "Number of acquired records = " << uRecords << endl;

  U32 TriggerAddress;
  U32 TimeStampHighPart;
  U32 TimeStampLowPart;
  double sampsPerSec = (double)getSamplingRate();



  //
  // Now that the zero, one, or more triggers have been saved into
  // the PCI-board memory, we must transfer them, one-by-one, 
  // to computer memory and then to hard disk.
  //
  // Now count the number of "valid" triggers
  // and then transfer the data for each acquisition 
  bool trigValid=false;
  int nValidTriggers = 0;

  long recordPreSize = -(long)sc->board()->getRecordPreSize();
  U32 recordSize     = sc->board()->getRecordSize();

  // Parameters used to convert bin number to time
  float dtScope = 1.0/getSamplingRate();
  int nbins     = int(getScopeConfig()->board()->getRecordSize());
  int nbinsPre  = int(getScopeConfig()->board()->getRecordPreSize());
  float xlow    = -nbinsPre*dtScope;
  float xup     = (nbins-nbinsPre)*dtScope;

  // Parameters used to convert ADC counts to volts
  float levelsToVoltsChA = getVoltageStep(chAId);
  float levelsToVoltsChB = getVoltageStep(chBId);
  float vZero = float(levelOfVolts(chAId, 0.0)); 
  //float vZero=128.0f; // valid for the ATS860, an 8-bit device (256 levels).

  //cout << "in ScopeAlazarATS860" << endl;
  //cout << "dtScope = " << dtScope << endl;
  //cout << "nbins, nbinsPre = " << nbins << ", " << nbinsPre << endl;
  //cout << "xlow, xup = " << xlow << ", " << xup << endl;

  U32 triggerCause;
  for (U32 recordNumber=1; recordNumber <= uRecords; recordNumber++) {

    triggerCause = AlazarGetWhoTriggeredBySystemHandle(hmaster, 1, recordNumber);
    //cout << "recordNumber, triggerCause = " << recordNumber << ", " << triggerCause << endl;
    // if triggered by CHANNEL_A or CHANNEL_B
    // (there are many more options that would be valid)
    // 1 = Channel A
    // 2 = Channel B
    // 4 = A and B
    trigValid = false;
    if (triggerCause==1 || triggerCause==2 || triggerCause==4) trigValid = true;
    //cout << "triggerCause, trigValid = " << triggerCause << ", " << trigValid << endl;
    
    if (!trigValid) continue;


    //cout << "reading record " << recordNumber << " from handle = " << hh << ", " << endl;

    nValidTriggers++;
    AlazarRead(hh, CHANNEL_A, _buffers[chAId], 
	       elementSize, recordNumber,
	       recordPreSize, recordSize); 
    
    status = AlazarGetTriggerAddress(hmaster, recordNumber, 
				     &TriggerAddress, 
				     &TimeStampHighPart,
				     &TimeStampLowPart);

    assert(apiError(status)==0);
    double triggerTime = timestampToSeconds(TimeStampHighPart,
					    TimeStampLowPart,
					    sampsPerSec);

    dataChanA->addTime( triggerTime );
    TString histoNameA = "scope_";
    histoNameA += _boardID;
    histoNameA += "_A_";
    histoNameA += nValidTriggers;

    // really want to save w/ x-axis = time and y-axis = voltage
    wfAId = dataChanA->addWaveform( new ScopeWaveformData(histoNameA, (int)sc->board()->getRecordSize(), xlow, xup,triggerTime,levelsToVoltsChA,vZero) );
    for (int ii=1; ii<=nbins; ii++) {
      dataChanA->wf(wfAId)->SetBinContent(ii, char(*(_buffers[chAId] + (ii-1) )));
    }

    AlazarRead(hh, CHANNEL_B, _buffers[chBId], 
	       elementSize, recordNumber,
	       recordPreSize, recordSize); 

    TString histoNameB = "scope_";
    histoNameB += _boardID;
    histoNameB += "_B_";
    histoNameB += nValidTriggers;
    
    wfBId = dataChanB->addWaveform( new ScopeWaveformData(histoNameB, (int)sc->board()->getRecordSize(), xlow, xup, triggerTime, levelsToVoltsChB,vZero) );
    for (int ii=1; ii<=nbins; ii++) {
      dataChanB->wf(wfBId)->SetBinContent(ii, char(*(_buffers[chBId] + (ii-1) )));
    }
    //cout << dataChanB->wf(wfBId)->GetName() << "   " << dataChanB->wf(wfBId)->GetMean() << endl;

  }

  data->setNValidTriggers(nValidTriggers);

  return 0;
}

int
ScopeAlazarATS860::configureHardwareForDAQ(ScopeConfig* sc) {

  //cout << GetName() << ".configureHardwareForDAQ()" << endl;

  setRecordSize(sc);
  setRecordCount(sc);
  setCaptureClock(sc);
  setChannel(CHANNEL_A, sc);  
  setChannel(CHANNEL_B, sc);
  setTrigger(sc);

  return 0;
}

// Deprecated... see setTrigLevel()
void 
ScopeAlazarATS860::setTriggerLevel1(U32 level) {
  // update the scope config
  ScopeConfig* sc = getScopeConfig();
  sc->trig()->setTrigLevel1(level);
  // update the hardware
  setTrigger(sc);
}
// Deprecated... see getTrigLevel()
U32
ScopeAlazarATS860::getTriggerLevel1() {
  cout << GetName() << ":  scopeConfig = " << getScopeConfig() << endl;

  return getScopeConfig()->trig()->trigLevel1();
}


U32
ScopeAlazarATS860::getTriggerLevel2() { 
  cout << GetName() << ":  scopeConfig = " << getScopeConfig() << endl;

  return getScopeConfig()->trig()->trigLevel2(); 
}

void
ScopeAlazarATS860::setTriggerLevel2(U32 level) {
  // update the scope config
  ScopeConfig* sc = getScopeConfig();
  sc->trig()->setTrigLevel2(level);
  // update the hardware
  setTrigger(sc);
}


void
ScopeAlazarATS860::setTriggerLevel1Volts(int chanNumber, float level_volts) {
  // convert volts to a digital level (0-255)
  U32 level = levelOfVolts(chanNumber, level_volts);
  // send that digital level to setTriggerLevel1
  setTriggerLevel1(level);

  cout << "set trigger to [volts] = " << voltsOfLevel(chanNumber, level) << endl;

}
void
ScopeAlazarATS860::setTriggerLevel2Volts(int chanNumber, float level_volts) {
  // convert volts to a digital level (0-255)
  U32 level = levelOfVolts(chanNumber, level_volts);
  // send that digital level to setTriggerLevel1
  setTriggerLevel2(level);

  cout << "set trigger to [volts] = " << voltsOfLevel(chanNumber, level) << endl;

}

U32
ScopeAlazarATS860::levelOfVolts(int chanNumber, float level_volts) {

  int level = int((level_volts-getVoltageMin(chanNumber))/getVoltageStep(chanNumber));

  if (level<0) level=0;
  if (level > getNLevels()) level=getNLevels()-1;
    
  return U32(level);
}

float
ScopeAlazarATS860::voltsOfLevel(int chanNumber, U32 level) {

  float level_volts = getVoltageMin(chanNumber) + getVoltageStep(chanNumber)*level;
  return level_volts;
}



void
ScopeAlazarATS860::setVoltageRange(int chanNumber, U32 range) { 
  // update the scope config
  ScopeConfig* sc = getScopeConfig();
  sc->chan(chanNumber)->setVoltageRange(range);
  // update the hardware
  U8 chanId;
  if (chanNumber == 0) 
    chanId = CHANNEL_A;
  else if (chanNumber == 1)
    chanId = CHANNEL_B;
  else {
    //cout << GetName() << ".setVoltageRange() ERROR\n" 
    // << "Invalid channel number provided:  chanNumber = " 
    // << chanNumber << endl;
    assert(!"Invalid channel provided to setVoltageRange()");
  }
  setChannel(chanId, sc);  
}

U32
ScopeAlazarATS860::getVoltageRange(int chanNumber) { return getScopeConfig()->chan(chanNumber)->getVoltageRange(); }

int
ScopeAlazarATS860::setTrigger(ScopeConfig* sc) {
  //cout << GetName() << ".setupTrigger()" << endl;

  status = AlazarSetTriggerOperation( getScopeHandle(),
				      sc->trig()->trigOp(),
				      sc->trig()->trigEn1(),
				      sc->trig()->trigSrc1(),
				      sc->trig()->trigSlope1(),
				      sc->trig()->trigLevel1(),
				      sc->trig()->trigEn2(),
				      sc->trig()->trigSrc2(),
				      sc->trig()->trigSlope2(),
				      sc->trig()->trigLevel2());
  assert(apiError(status)==0);

  // set up external trigger if necessary
  // Note, there is no API error checking on this yet...
  if ( (sc->trig()->trigSrc1() == TRIG_EXTERNAL) ||
       (sc->trig()->trigSrc2() == TRIG_EXTERNAL) ) {
    status = AlazarSetExternalTrigger(getScopeHandle(), DC_COUPLING, ETR_X1);
    assert(apiError(status)==0);
  }

  // set up a trigger timeout
  //U32 timeout_us = 1000000;
  //status = AlazarSetTriggerTimeOut(getScopeHandle(), timeout_us);
  //assert(apiError(status)==0);

  //return apiError(status);
  return 0;
}

U32
ScopeAlazarATS860::setTriggerDelay(U32 delay) {
  
  return AlazarSetTriggerDelay(getScopeHandle(), delay);
}

U32
ScopeAlazarATS860::setTriggerSlope(U32 slope) {
  cout << GetName() << "::setTriggerSlope not done" << endl;
  return 0;
}


/////////////////////////////////////////////////////////////////
// Start of trigger set/get methods
void
ScopeAlazarATS860::setTrigEngineOperation(ScopeTriggerEngOp operation) {
  // Specify the logic criteria between the two trigger engines
  // operation:  the logic criteria
  //  valid values are:
  //    TRIG_ENGINE_OP_J
  //    TRIG_ENGINE_OP_K
  //    TRIG_ENGINE_OP_J_OR_K
  //    TRIG_ENGINE_OP_J_AND_K
  //    TRIG_ENGINE_OP_J_XOR_K
  //    TRIG_ENGINE_OP_J_AND_NOT_K
  //    TRIG_ENGINE_OP_NOT_J_AND_K

  // update the scope config
  ScopeConfig* sc = getScopeConfig();
  sc->trig()->setTrigOp(GenericToAlazarEngOp(operation));
  
  // update the hardware
  setTrigger(sc);
}
ScopeTriggerEngOp ScopeAlazarATS860::getTrigEngineOperation() { return AlazarToGenericEngOp(getScopeConfig()->trig()->trigOp()); }

void
ScopeAlazarATS860::setTrigLevel(int itrg, U32 level) {
  // itrg:   1 or 2 (which trigger engine are you configuring)
  // level:  The trigger level value (0..255) for trigger engine itrg.
  //         128 is zero percent, 255 is +100%, 0 is -100%
  
  // update the scope config
  ScopeConfig* sc = getScopeConfig();
  if (itrg == 1) sc->trig()->setTrigLevel1(level);
  else if (itrg == 2) sc->trig()->setTrigLevel2(level);
  else assert(!"setTrigLevel: Invalid trigger number, itrg");
  // update the hardware
  setTrigger(sc);
}
U32 ScopeAlazarATS860::getTrigLevel(int itrg) { 
  if (itrg == 1) return getScopeConfig()->trig()->trigLevel1(); 
  else if (itrg == 2) return getScopeConfig()->trig()->trigLevel2(); 
  else assert(!"getTrigLevel: Invalid trigger number, itrg");
}

void 
ScopeAlazarATS860::setTrigEngine(int itrg, ScopeTriggerEngine engine) {
  // Set up the trigger engine
  // itrg:   1 or 2 (which trigger engine are you configuring)
  // engine options are:
  //  TRIG_ENGINE_J
  //  TRIG_ENGINE_K

  // update the scope config
  ScopeConfig* sc = getScopeConfig();

  AU32 engineAlazar = GenericToAlazarEngine(engine);

  if (itrg == 1) sc->trig()->setTrigEngine1(engineAlazar);
  else if (itrg == 2) sc->trig()->setTrigEngine2(engineAlazar);
  else assert(!"setTrigEngine: Invalid trigger number, itrg");
  
  // update the hardware
  setTrigger(sc);
}
ScopeTriggerEngine ScopeAlazarATS860::getTrigEngine(int itrg) {  
  AU32 engineAlazar;
  if (itrg == 1) engineAlazar = getScopeConfig()->trig()->trigEn1(); 
  else if (itrg == 2) engineAlazar = getScopeConfig()->trig()->trigEn2(); 
  else assert(!"getTrigEngine: Invalid trigger number, itrg");

  return AlazarToGenericEngine(engineAlazar);
}
  
void 
ScopeAlazarATS860::setTrigSource(int itrg, ScopeTriggerSource source) {
  // itrg:   1 or 2 (which trigger engine are you configuring)
  // source:  specifies which source to use as trigger
  //  valid values are:  
  //    TRIG_CHAN_A
  //    TRIG_CHAN_B
  //    TRIG_EXTERNAL
  //    TRIG_DISABLE

  // update the scope config
  ScopeConfig* sc = getScopeConfig();

  U32 sourceAlazar = GenericToAlazarSource(source);

  if (itrg == 1) sc->trig()->setTrigSource1(sourceAlazar);
  else if (itrg == 2) sc->trig()->setTrigSource2(sourceAlazar);
  else assert(!"setTrigSource: Invalid trigger number, itrg");

  // update the hardware
  setTrigger(sc);
}
ScopeTriggerSource ScopeAlazarATS860::getTrigSource(int itrg) { 
  AU32 sourceAlazar;
  if (itrg == 1) sourceAlazar = getScopeConfig()->trig()->trigSrc1(); 
  else if (itrg == 2) sourceAlazar = getScopeConfig()->trig()->trigSrc2(); 
  else assert(!"getTrigEngine: Invalid trigger number, itrg");

  return AlazarToGenericSource(sourceAlazar);

}

/////////////////////////////////////////////////////////////////////////
// type conversion utilities
ScopeTriggerSource 
ScopeAlazarATS860::AlazarToGenericSource(AU32 source) {
  // convert between the Alazar type for Source and the ScopeType.hh value
  if (source == TRIG_DISABLE) return SCOPE_TRIGGER_SOURCE_DISABLE;
  else if (source == TRIG_EXTERNAL) return SCOPE_TRIGGER_SOURCE_EXTERNAL;
  else if (source == TRIG_CHAN_A) return SCOPE_TRIGGER_SOURCE_CHAN_A;
  else if (source == TRIG_CHAN_B) return SCOPE_TRIGGER_SOURCE_CHAN_B;
  else assert(!"AlazarToGenericSource: invalid AU32 source value");
}
AU32 
ScopeAlazarATS860::GenericToAlazarSource(ScopeTriggerSource source) {
  if (source == SCOPE_TRIGGER_SOURCE_DISABLE) return TRIG_DISABLE;
  else if (source == SCOPE_TRIGGER_SOURCE_EXTERNAL) return TRIG_EXTERNAL;
  else if (source == SCOPE_TRIGGER_SOURCE_CHAN_A) return TRIG_CHAN_A;
  else if (source == SCOPE_TRIGGER_SOURCE_CHAN_B) return TRIG_CHAN_B;
  else assert(!"GenericToAlazarSource: invalid ScopeTriggerSource value");
}
ScopeTriggerSlope
ScopeAlazarATS860::AlazarToGenericSlope(AU32 slope) {
  if (slope == TRIGGER_SLOPE_POSITIVE) return SCOPE_TRIGGER_SLOPE_POSITIVE;
  else if (slope == TRIGGER_SLOPE_NEGATIVE) return SCOPE_TRIGGER_SLOPE_NEGATIVE;
  else assert(!"AlazarToGenericSlope: Invalid input AU32 slope value");
}
AU32
ScopeAlazarATS860::GenericToAlazarSlope(ScopeTriggerSlope slope) {
  if (slope == SCOPE_TRIGGER_SLOPE_POSITIVE) return TRIGGER_SLOPE_POSITIVE;
  else if (slope == SCOPE_TRIGGER_SLOPE_NEGATIVE) return TRIGGER_SLOPE_NEGATIVE;
  else assert(!"GenericToAlazarSlope: Invalid input ScopeTriggerSlope slope value");
}
ScopeTriggerEngine
ScopeAlazarATS860::AlazarToGenericEngine(AU32 engine) {
  if (engine == TRIG_ENGINE_J) return SCOPE_TRIGGER_ENGINE_J;
  else if (engine == TRIG_ENGINE_K) return SCOPE_TRIGGER_ENGINE_K;
  else assert(!"AlazarToGenericEngine:  Invalid input AU32 engine");
}
AU32
ScopeAlazarATS860::GenericToAlazarEngine(ScopeTriggerEngine engine) {
  if (engine == SCOPE_TRIGGER_ENGINE_J) return TRIG_ENGINE_J;
  else if (engine == SCOPE_TRIGGER_ENGINE_K) return TRIG_ENGINE_K;
  else assert(!"GenericToAlazarEngine:  Invalid input ScopeTriggerEngine engine");
}

ScopeTriggerEngOp
ScopeAlazarATS860::AlazarToGenericEngOp(AU32 engOp) {
  if (engOp == TRIG_ENGINE_OP_J) return SCOPE_TRIGGER_ENG_OP_J;
  else if (engOp == TRIG_ENGINE_OP_K) return SCOPE_TRIGGER_ENG_OP_K;
  else if (engOp == TRIG_ENGINE_OP_J_OR_K) return SCOPE_TRIGGER_ENG_OP_J_OR_K;
  else if (engOp == TRIG_ENGINE_OP_J_AND_K) return SCOPE_TRIGGER_ENG_OP_J_AND_K;
  else if (engOp == TRIG_ENGINE_OP_J_XOR_K) return SCOPE_TRIGGER_ENG_OP_J_XOR_K;
  else if (engOp == TRIG_ENGINE_OP_J_AND_NOT_K) return SCOPE_TRIGGER_ENG_OP_J_AND_NOT_K;
  else if (engOp == TRIG_ENGINE_OP_NOT_J_AND_K) return SCOPE_TRIGGER_ENG_OP_NOT_J_AND_K;
  else assert(!"AlazarToGenericEngOp:  Invalid input AU32 engOp");
}
AU32
ScopeAlazarATS860::GenericToAlazarEngOp(ScopeTriggerEngOp engOp) {
  if (engOp == SCOPE_TRIGGER_ENG_OP_J) return TRIG_ENGINE_OP_J;
  else if (engOp == SCOPE_TRIGGER_ENG_OP_K) return TRIG_ENGINE_OP_K;
  else if (engOp == SCOPE_TRIGGER_ENG_OP_J_OR_K) return TRIG_ENGINE_OP_J_OR_K;
  else if (engOp == SCOPE_TRIGGER_ENG_OP_J_AND_K) return TRIG_ENGINE_OP_J_AND_K;
  else if (engOp == SCOPE_TRIGGER_ENG_OP_J_XOR_K) return TRIG_ENGINE_OP_J_XOR_K;
  else if (engOp == SCOPE_TRIGGER_ENG_OP_J_AND_NOT_K) return TRIG_ENGINE_OP_J_AND_NOT_K;
  else if (engOp == SCOPE_TRIGGER_ENG_OP_NOT_J_AND_K) return TRIG_ENGINE_OP_NOT_J_AND_K;
  else assert(!"GenericToAlazarEngOp:  Invalid input ScopeTriggerEngOp engOp");
}

// end of type conversion utilities
/////////////////////////////////////////////////////////////////////////

  
void
ScopeAlazarATS860::setTrigSlope(int itrg, ScopeTriggerSlope slope) {
  // itrg:   1 or 2 (which trigger engine are you configuring)
  // slope:  specify the slope criteria for the trigger
  //  valid values are: 
  //    SCOPE_TRIGGER_SLOPE_POSITIVE
  //    SCOPE_TRIGGER_SLOPE_NEGATIVE

  // update the scope config
  ScopeConfig* sc = getScopeConfig();
  U32 slopeAlazar = GenericToAlazarSlope(slope);

  if (itrg == 1) sc->trig()->setTrigSlope1(slopeAlazar);
  else if (itrg == 2) sc->trig()->setTrigSlope2(slopeAlazar);
  else assert(!"setTrigSlope: Invalid trigger number, itrg");
  
  // update the hardware
  setTrigger(sc);
}
ScopeTriggerSlope ScopeAlazarATS860::getTrigSlope(int itrg) { 
  U32 slopeAlazar;
  if (itrg == 1) slopeAlazar = getScopeConfig()->trig()->trigSlope1(); 
  else if (itrg == 2) slopeAlazar = getScopeConfig()->trig()->trigSlope2(); 
  else assert(!"getTrigSlope: Invalid trigger number, itrg");

  return AlazarToGenericSlope(slopeAlazar);
}
  
int
ScopeAlazarATS860::setChannel(U8 chanId, ScopeConfig* sc) {
  //cout << GetName() << ".setChannel() " << (int)chanId << endl;
  //cout << "     CHANNEL_A = " << (int)CHANNEL_A << endl;
  //cout << "     CHANNEL_B = " << (int)CHANNEL_B << endl;
  // chanId is the channel label specified by Alazar
  // e.g. CHANNEL_A or CHANNEL_B

  // "channel" is an index into the vector of channels in ScopeConfig
  U8 channel;

  // this relies on addChannelConfig() being added in the order: 
  //  CHANNEL_A, CHANNEL_B
  // Should add capability in ScopeConfig class to get channel id by
  // label.  Then you can say: 
  // channel = sc.chanIdOf(CHANNEL_A)
  if (chanId == CHANNEL_A)
    channel = 0;
  else if (chanId == CHANNEL_B)
    channel = 1;
  else {
    cout << GetName() << ".setChannel() ERROR\n" 
	 << "Invalid channel id provided:  chanId = " << chanId << endl;
    assert(!"Invalid channel provided to setChannel()");
  }

  //cout << "getScopeHandle() = " << getScopeHandle() << endl;
  //cout << "channel = " << int(channel) << endl;
  //cout << "sc->chan(" << int(channel) << ")->getCoupling() = " << sc->chan(channel)->getCoupling() << endl;
  //cout << "sc->chan(channel)->getVoltageRange() = " << sc->chan(channel)->getVoltageRange() << endl;
  //cout << "sc->chan(channel)->getInputImpedance() = " << sc->chan(channel)->getInputImpedance() << endl;
  //cout << "status = " << status << endl;
  status = AlazarInputControl( getScopeHandle(), chanId,
  			       sc->chan(channel)->getCoupling(),
  			       sc->chan(channel)->getVoltageRange(),
  			       sc->chan(channel)->getInputImpedance());
  //cout << "status = " << status << endl;
  assert(apiError(status)==0);
  //cout << "getInputImpedance() = " << getInputImpedance(channel) << endl;

  return 0;
}

int
ScopeAlazarATS860::setCaptureClock(ScopeConfig* sc) {
  //cout << GetName() << ".setCaptureClock()" << endl;
  status = AlazarSetCaptureClock( getScopeHandle(), 
				  sc->clock()->getClockSource(),
				  sc->clock()->getClockRate(),
				  sc->clock()->getClockEdge(),
				  0);
  assert(apiError(status)==0);
  return 0;
}

int
ScopeAlazarATS860::setRecordSize(ScopeConfig* sc) {
  //cout << GetName() << ".setRecordSize()" << endl;
  status = AlazarSetRecordSize( getScopeHandle(), 
				sc->board()->getRecordPreSize(),
				sc->board()->getRecordPostSize() );
  assert(apiError(status)==0);
  return 0;
}

int
ScopeAlazarATS860::setRecordCount(ScopeConfig* sc) {
  //cout << GetName() << ".setRecordCount()" << endl;
  status = AlazarSetRecordCount( getScopeHandle(), 
				 sc->board()->getNRecordsPerBuffer() );
  assert(apiError(status)==0);
  return 0;
}


ScopeConfig*
ScopeAlazarATS860::makeScopeConfig(BoardDef* bd, 
				   TString* chALabel, 
				   TString* chBLabel) {
  // Set up information about the board
  ScopeBoardInfo* boardInfo;
  boardInfo = new ScopeBoardInfo();
  boardInfo->setBoardType(&Scope::ALAZAR_ATS860);
  boardInfo->setBoardChannelNumbers(ALAZAR_ATS860_CHANNELS, 2);
  //cout << "bd->PreDepth  = " << (int)bd->PreDepth << endl;
  //cout << "bd->RecLength = " << (int)bd->RecLength << endl;
  boardInfo->setRecordPreSize(bd->PreDepth);
  boardInfo->setRecordPostSize(bd->RecLength-bd->PreDepth);
  //cout << "boardInfo->getRecordSize() = " 
  //     << (int)boardInfo->getRecordSize() << endl;
  boardInfo->setBytesPerRecord(computeBytesPerRecord(bd));
  boardInfo->setNRecordsPerBuffer(bd->RecordCount);
  //boardInfo->setNRecordsPerBuffer(N_RECORDS_PER_BUFFER);

  // Set up the clock
  ScopeClockConfig* clockConfig;
  clockConfig = new ScopeClockConfig(bd->ClockSource, bd->SampleRate,
                                     bd->ClockEdge, 0);
  // Set up the trigger
  ScopeTriggerConfig* triggerConfig;
  triggerConfig = new ScopeTriggerConfig(bd->TriEngOperation,
                                         bd->TriggerEngine1, 
					 bd->TrigEngSource1,
                                         bd->TrigEngSlope1,  
					 bd->TrigEngLevel1,
                                         bd->TriggerEngine2, 
					 bd->TrigEngSource2,
                                         bd->TrigEngSlope2,  
					 bd->TrigEngLevel2);
  //Set up each channel
  ScopeChannelConfig* chAConfig;
  ScopeChannelConfig* chBConfig;
  chAConfig = new ScopeChannelConfig(CHANNEL_A, bd->CouplingChanA,
                                     bd->InputRangeChanA,
                                     bd->InputImpedChanA, chALabel);
  chBConfig = new ScopeChannelConfig(CHANNEL_B, bd->CouplingChanB,
                                     bd->InputRangeChanB,
                                     bd->InputImpedChanB, chBLabel);
  // Make and populate the master config object
  ScopeConfig* scopeConfig = new ScopeConfig();
  //cout << "bd->PreDepth  = " << (int)bd->PreDepth << endl;
  //cout << "bd->RecLength = " << (int)bd->RecLength << endl;
  //scopeConfig->setRecordPreSize(bd->PreDepth);
  //scopeConfig->setRecordPostSize(bd->RecLength-bd->PreDepth);
  //cout << "scopeConfig->getRecordSize() = " 
  //<< (int)scopeConfig->getRecordSize() << endl;
  //scopeConfig->setBytesPerRecord(computeBytesPerRecord(bd));
  //scopeConfig->setNRecordsPerBuffer(N_RECORDS_PER_BUFFER);
  //scopeConfig->setNRecordsPerBuffer(bd->RecordCount);
  scopeConfig->setBoardConfig(boardInfo);
  scopeConfig->setClockConfig(clockConfig);
  scopeConfig->setTriggerConfig(triggerConfig);
  scopeConfig->addChannelConfig(chAConfig);
  scopeConfig->addChannelConfig(chBConfig);

  //cout << "scopeConfig->getRecordSize() = " << (int)scopeConfig->board()->getRecordSize() << endl;

  return scopeConfig;
}

////  DEFUNCT -- Do not use 
//float
//ScopeAlazarATS860::getDtOfClockRate(AU32 sampleRate) {
//  float dtPerSample;
//
//  if (sampleRate == SAMPLE_RATE_1KSPS) {
//    dtPerSample = 1e-3;
//  } else if (sampleRate == SAMPLE_RATE_2KSPS) {
//    dtPerSample = 5e-4;
//  } else if (sampleRate == SAMPLE_RATE_5KSPS) {
//    dtPerSample = 2e-4;
//  } else if (sampleRate == SAMPLE_RATE_10KSPS) {
//    dtPerSample = 1e-4;
//  } else if (sampleRate == SAMPLE_RATE_20KSPS) {
//    dtPerSample = 5e-5;
//  } else if (sampleRate == SAMPLE_RATE_50KSPS) {
//    dtPerSample = 2e-5;
//  } else if (sampleRate == SAMPLE_RATE_100KSPS) {
//    dtPerSample = 1e-5;
//  } else if (sampleRate == SAMPLE_RATE_200KSPS) {
//    dtPerSample = 5e-6;
//  } else if (sampleRate == SAMPLE_RATE_500KSPS) {
//    dtPerSample = 2e-6;
//  } else if (sampleRate == SAMPLE_RATE_1MSPS) {
//    dtPerSample = 1e-6;
//  } else if (sampleRate == SAMPLE_RATE_2MSPS) {
//    dtPerSample = 5e-7;
//  } else if (sampleRate == SAMPLE_RATE_5MSPS) {
//    dtPerSample = 2e-7;
//  } else if (sampleRate == SAMPLE_RATE_10MSPS) {
//    dtPerSample = 1e-7;
//  } else if (sampleRate == SAMPLE_RATE_20MSPS) {
//    dtPerSample = 5e-8;
//  } else if (sampleRate == SAMPLE_RATE_25MSPS) {
//    dtPerSample = 4e-8;
//  } else if (sampleRate == SAMPLE_RATE_50MSPS) {
//    dtPerSample = 2e-8;
//  } else if (sampleRate == SAMPLE_RATE_100MSPS) {
//    dtPerSample = 1e-8;
//  } else if (sampleRate == SAMPLE_RATE_125MSPS) {
//    dtPerSample = 8e-9;
//  } else if (sampleRate == SAMPLE_RATE_200MSPS) {
//    dtPerSample = 5e-9;
//  } else if (sampleRate == SAMPLE_RATE_250MSPS) {
//    dtPerSample = 4e-9;
//  } else if (sampleRate == SAMPLE_RATE_500MSPS) {
//    dtPerSample = 2e-9;
//  } else if (sampleRate == SAMPLE_RATE_1GSPS) {
//    dtPerSample = 1e-9;
//  } else {
//    cout << GetName() 
//	 << ".getDtOfClockRate(): unrecognized sampleRate" << endl;
//    assert(!"error: sampleRate not recognized");
//  }
//
//  return dtPerSample;
//}

float
ScopeAlazarATS860::getSamplingRate() {
  //ScopeConfig* sc = getScopeConfig();
  //AU32 sampleRate = sc->clock()->getClockRate();

  AU32 sampleRate = getScopeConfig()->clock()->getClockRate();

  float samplingRate;

  if (sampleRate == SAMPLE_RATE_1KSPS) {
    samplingRate = 1e3;
  } else if (sampleRate == SAMPLE_RATE_2KSPS) {
    samplingRate = 2e3;
  } else if (sampleRate == SAMPLE_RATE_5KSPS) {
    samplingRate = 5e3;
  } else if (sampleRate == SAMPLE_RATE_10KSPS) {
    samplingRate = 10e3;
  } else if (sampleRate == SAMPLE_RATE_20KSPS) {
    samplingRate = 20e3;
  } else if (sampleRate == SAMPLE_RATE_50KSPS) {
    samplingRate = 50e3;
  } else if (sampleRate == SAMPLE_RATE_100KSPS) {
    samplingRate = 100e3;
  } else if (sampleRate == SAMPLE_RATE_200KSPS) {
    samplingRate = 200e3;
  } else if (sampleRate == SAMPLE_RATE_500KSPS) {
    samplingRate = 500e3;
  } else if (sampleRate == SAMPLE_RATE_1MSPS) {
    samplingRate = 1e6;
  } else if (sampleRate == SAMPLE_RATE_2MSPS) {
    samplingRate = 2e6;
  } else if (sampleRate == SAMPLE_RATE_5MSPS) {
    samplingRate = 5e6;
  } else if (sampleRate == SAMPLE_RATE_10MSPS) {
    samplingRate = 10e6;
  } else if (sampleRate == SAMPLE_RATE_20MSPS) {
    samplingRate = 20e6;
  } else if (sampleRate == SAMPLE_RATE_25MSPS) {
    samplingRate = 25e6;
  } else if (sampleRate == SAMPLE_RATE_50MSPS) {
    samplingRate = 50e6;
  } else if (sampleRate == SAMPLE_RATE_100MSPS) {
    samplingRate = 100e6;
  } else if (sampleRate == SAMPLE_RATE_125MSPS) {
    samplingRate = 125e6;
  } else if (sampleRate == SAMPLE_RATE_200MSPS) {
    samplingRate = 200e6;
  } else if (sampleRate == SAMPLE_RATE_250MSPS) {
    samplingRate = 250e6;
  } else if (sampleRate == SAMPLE_RATE_500MSPS) {
    samplingRate = 500e6;
  } else if (sampleRate == SAMPLE_RATE_1GSPS) {
    samplingRate = 1e9;
  } else {
    cout << GetName() 
	 << ".getSamplingRate() unrecognized sampleRate" << endl;
    assert(!"error: sampleRate not recognized");
  }

  return samplingRate;
}

float
ScopeAlazarATS860::getVoltageMax(int chanNum) {
  
  AU32 voltRangeSpecifier = getScopeConfig()->chan(chanNum)->getVoltageRange();

  if (voltRangeSpecifier == INPUT_RANGE_PM_20_MV) {
    return 0.02f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_40_MV) {
    return 0.04f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_50_MV) {
    return 0.05f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_80_MV) {
    return 0.08f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_100_MV) {
    return 0.1f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_200_MV) {
    return 0.2f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_400_MV) {
    return 0.4f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_500_MV) {
    return 0.5f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_800_MV) {
    return 0.8f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_1_V) {
    return 1.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_2_V) {
    return 2.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_4_V) {
    return 4.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_5_V) {
    return 5.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_8_V) {
    return 8.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_10_V) {
    return 10.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_20_V) {
    return 20.0f;
  } else {
    cout << GetName() << ".getVoltageMax(): unrecognized voltRangeSpecifier" << endl;
    assert(!"error: voltRangeSpecifier not recognized");
  }
}


float
ScopeAlazarATS860::getVoltageMin(int chanNum) {
  return -getVoltageMax(chanNum);
}


void
ScopeAlazarATS860::setInputCoupling(int chanNumber, ScopeInputCoupling coupling) { 
  // chanNumber = 0 (CHANNEL_A)
  //              1 (CHANNEL_B)
  // coupling   = 1 (SCOPE_COUPLING_AC)  (see ScopeTypes.hh)
  //              2 (SCOPE_COUPLING_DC)

  // update the scope config
  ScopeConfig* sc = getScopeConfig();

  U32 alazarCoupling;
  if (coupling == SCOPE_COUPLING_AC) alazarCoupling = AC_COUPLING;
  else if (coupling == SCOPE_COUPLING_DC) alazarCoupling = DC_COUPLING;
  else {
    assert(!"Invalid coupling provided to setInputCoupling()");
  }
  sc->chan(chanNumber)->setInputCoupling(alazarCoupling);
  // update the hardware
  U8 chanId;
  if (chanNumber == 0) 
    chanId = CHANNEL_A;
  else if (chanNumber == 1)
    chanId = CHANNEL_B;
  else {
    assert(!"Invalid channel provided to setInputCoupling()");
  }
  setChannel(chanId, sc);  
}

ScopeInputCoupling
ScopeAlazarATS860::getInputCoupling(int chanNum) {
  AU32 inputCoupling = getScopeConfig()->chan(chanNum)->getCoupling();
  if (inputCoupling == AC_COUPLING) {
    return SCOPE_COUPLING_AC;
  } else if (inputCoupling == DC_COUPLING) {
    return SCOPE_COUPLING_DC;
  } else {
    cout << GetName() << ".getInputCoupling(): unrecognized inputCoupling" << endl;
    assert(!"error: inputCoupling not recognized");
  }
}


void
ScopeAlazarATS860::setInputImpedance(int chanNumber, U32 impedance) { 
  // chanNumber = 0 (CHANNEL_A)
  //              1 (CHANNEL_B)
  // impedance  = 1 (IMPEDANCE_1M_OHM)
  //              2 (IMPEDANCE_50_OHM)

  // update the scope config
  ScopeConfig* sc = getScopeConfig();
  sc->chan(chanNumber)->setInputImpedance(impedance);
  // update the hardware
  U8 chanId;
  if (chanNumber == 0) 
    chanId = CHANNEL_A;
  else if (chanNumber == 1)
    chanId = CHANNEL_B;
  else {
    assert(!"Invalid channel provided to setInputImpedance()");
  }
  setChannel(chanId, sc);  
}

float
ScopeAlazarATS860::getInputImpedance(int chanNum) {
  AU32 inputImpedance = getScopeConfig()->chan(chanNum)->getInputImpedance();
  if (inputImpedance == IMPEDANCE_1M_OHM) {
    return 1e6;
  } else if (inputImpedance == IMPEDANCE_50_OHM) {
    return 50.0;
  } else {
    cout << GetName() << ".getInputImpedance(): unrecognized inputImpedance" << endl;
    assert(!"error: inputImpedance not recognized");
  }
}

float ScopeAlazarATS860::getVoltageStep(int chanNum) { 
  int nLevels = getNLevels();

  // could replace with getVoltageRange()
  float vMin = getVoltageMin(chanNum);
  float vMax = getVoltageMax(chanNum);
  
  return (vMax-vMin)/nLevels;
}

float ScopeAlazarATS860::getTriggerLevel(int chanNum) {
  // For the ATS860, a trigger level ranges from 0 to 255
  // 128 = 0% 
  // 255 = +100% of full range
  //   0 = -100% of full range

  // this is hard to do correctly because the trigger level is 
  // associated with a trigger "engine", but the source of a 
  // trigger engine could be either channel

  // The following is, in general, wrong.
  // But in the case where the source of trigger engine 1 is CHANNEL_A
  // then it is right.
  int triggerLevel = int(getScopeConfig()->trig()->trigLevel1());

  // convert the level (in ADU) to a level in volts
  int triggerZeroLevel = 128;
  if (triggerLevel == triggerZeroLevel) {
    return 0.0; 
  } else if (triggerLevel > triggerZeroLevel) {
    float fracRangePerLevel = 1.0/127;
    return fracRangePerLevel*(triggerLevel-triggerZeroLevel)*getVoltageMax(chanNum);
  } else if ( (triggerLevel < triggerZeroLevel) && (triggerLevel >= 0) ) {
    float fracRangePerLevel = 1.0/128;
    return fracRangePerLevel*(triggerZeroLevel-triggerLevel)*getVoltageMin(chanNum);
  } else {
    cout << GetName() << ".getTriggerLevel(): unrecognized triggerLevel" << endl;
    assert(!"error: triggerLevel not recognized");
  }

}

ScopeTriggerSlope ScopeAlazarATS860::getTriggerSlope(int chanNum) {
  // ******WARNING****** 
  // same "this might be wrong" caveat as getTriggerLevel()
  //

  AU32 triggerSlope = getScopeConfig()->trig()->trigSlope1();
  if (triggerSlope == TRIGGER_SLOPE_POSITIVE) {
    return SCOPE_TRIGGER_SLOPE_POSITIVE;
  } else if (triggerSlope == TRIGGER_SLOPE_NEGATIVE) {
    return SCOPE_TRIGGER_SLOPE_NEGATIVE;
  } else {
    cout << GetName() << ".getTriggerSlope(): unrecognized triggerSlope" << endl;
    assert(!"error: triggerSlope not recognized");
  }
}


// obsolete --- do not use -- see getVoltageMax instead
float
ScopeAlazarATS860::getMaxVoltage(U32 voltRangeSpecifier) {
  if (voltRangeSpecifier == INPUT_RANGE_PM_20_MV) {
    return 0.02f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_40_MV) {
    return 0.04f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_50_MV) {
    return 0.05f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_80_MV) {
    return 0.08f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_100_MV) {
    return 0.1f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_200_MV) {
    return 0.2f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_400_MV) {
    return 0.4f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_500_MV) {
    return 0.5f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_800_MV) {
    return 0.8f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_1_V) {
    return 1.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_2_V) {
    return 2.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_4_V) {
    return 4.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_5_V) {
    return 5.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_8_V) {
    return 8.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_10_V) {
    return 10.0f;
  } else if (voltRangeSpecifier == INPUT_RANGE_PM_20_V) {
    return 20.0f;
  } else {
    cout << GetName() << ".getMaxVoltage(): unrecognized voltRangeSpecifier" << endl;
    assert(!"error: voltRangeSpecifier not recognized");
  }
}

int
ScopeAlazarATS860::computeBytesPerRecord(BoardDef* bd) {
  int bytesPerRecord = 0;
  if (_acqMode == ACQ_MODE_TRADITIONAL)
    bytesPerRecord = (bd->RecLength + 16)*sizeof(U8);
  else if (_acqMode == ACQ_MODE_AUTODMA) {
    cout << GetName() << ".computeBytesPerRecord() ERROR\n"
	 << "AUTODMA mode not yet enabled" << endl;
    assert(!"ERROR: AUTODMA mode not yet enabled");
  }
  return bytesPerRecord;
}

BoardDef
ScopeAlazarATS860::getBoardDef(Scope::DAQ_State daqstate) {
  // Configure scope board for differrent readout schemes

  BoardDef bd;

  _state=daqstate;

  if (daqstate == Scope::DAQ_CHARGE) {
    bd.RecordCount     = N_RECORDS_PER_BUFFER;     // this could depend on _acqMode
    bd.RecLength       = 12*1024;
    bd.PreDepth        = 4*1024;
    bd.ClockSource     = INTERNAL_CLOCK;
    bd.ClockEdge       = CLOCK_EDGE_RISING;
    bd.SampleRate      = SAMPLE_RATE_250MSPS;  

    bd.CouplingChanA   = AC_COUPLING;
    bd.InputRangeChanA = INPUT_RANGE_PM_4_V; 
    bd.InputImpedChanA = IMPEDANCE_1M_OHM;

    bd.CouplingChanB   = AC_COUPLING;
    bd.InputRangeChanB = INPUT_RANGE_PM_4_V;
    bd.InputImpedChanB = IMPEDANCE_1M_OHM ;

    bd.TriEngOperation = TRIG_ENGINE_OP_J_OR_K; 

    bd.TriggerEngine1  = TRIG_ENGINE_J;
    bd.TrigEngSource1  = TRIG_CHAN_A;
    bd.TrigEngSlope1   = TRIGGER_SLOPE_POSITIVE;
    bd.TrigEngLevel1   = 128-10;                 

    bd.TriggerEngine2  = TRIG_ENGINE_K;
    bd.TrigEngSource2  = TRIG_CHAN_B;
    bd.TrigEngSlope2   = TRIGGER_SLOPE_POSITIVE;
    bd.TrigEngLevel2   = 128+10;

    if (_boardID) { // take triggers from 1st board only
      bd.TrigEngSource1 = TRIG_DISABLE;
      bd.TrigEngSource2 = TRIG_DISABLE;
    }

  } else if (daqstate==Scope::DAQ_TEST) {
    bd.RecordCount     = N_RECORDS_PER_BUFFER;     // this could depend on _acqMode
    bd.RecLength       = 4*1024;
    bd.PreDepth        = 1024;
    bd.ClockSource     = INTERNAL_CLOCK;
    bd.ClockEdge       = CLOCK_EDGE_RISING;
    bd.SampleRate      = SAMPLE_RATE_200KSPS;
    bd.SampleRate      = SAMPLE_RATE_1MSPS;
    //bd.SampleRate      = SAMPLE_RATE_5MSPS;
    bd.CouplingChanA   = DC_COUPLING;
    bd.InputRangeChanA = INPUT_RANGE_PM_10_V;
    bd.InputImpedChanA = IMPEDANCE_1M_OHM;
    bd.CouplingChanB   = bd.CouplingChanA;
    bd.InputRangeChanB = bd.InputRangeChanA;
    bd.InputImpedChanB = bd.InputImpedChanA;
    bd.TriEngOperation = TRIG_ENGINE_OP_J; // 0
    bd.TriggerEngine1  = TRIG_ENGINE_J;
    bd.TrigEngSource1  = TRIG_CHAN_A;
    bd.TrigEngSlope1   = TRIGGER_SLOPE_NEGATIVE;
    bd.TrigEngLevel1   = 128-15;  // negative going pulses
    //bd.TrigEngLevel1   = 128+15;  // positive going pulses
    //bd.TrigEngLevel1   = 128+127;  // positive going pulses
    bd.TriggerEngine2  = TRIG_ENGINE_K;
    bd.TrigEngSource2  = TRIG_DISABLE;
    bd.TrigEngSlope2   = TRIGGER_SLOPE_POSITIVE;
    bd.TrigEngLevel2   = 128;
  } else if (daqstate==Scope::DAQ_CV1) {
    bd.RecordCount     = N_RECORDS_PER_BUFFER;     // this could depend on _acqMode
    bd.RecLength       = 4*1024;
    bd.PreDepth        = 1*1024;
    bd.ClockSource     = INTERNAL_CLOCK;
    bd.ClockEdge       = CLOCK_EDGE_RISING;
    //bd.SampleRate      = SAMPLE_RATE_10MSPS;
    //bd.SampleRate      = SAMPLE_RATE_1MSPS;
    bd.SampleRate      = SAMPLE_RATE_5MSPS;
    bd.CouplingChanA   = DC_COUPLING;
    bd.InputRangeChanA = INPUT_RANGE_PM_100_MV; 
    //bd.InputRangeChanA = INPUT_RANGE_PM_10_V; 
    //bd.InputRangeChanA = INPUT_RANGE_PM_40_MV; // 90 of Alazar manual
    //bd.InputImpedChanA = IMPEDANCE_50_OHM;
    bd.InputImpedChanA = IMPEDANCE_1M_OHM;
    bd.CouplingChanB   = bd.CouplingChanA;
    bd.InputRangeChanB = bd.InputRangeChanA;
    bd.InputImpedChanB = bd.InputImpedChanA;
    bd.TriEngOperation = TRIG_ENGINE_OP_J; // 0
    bd.TriggerEngine1  = TRIG_ENGINE_J;
    bd.TrigEngSource1  = TRIG_CHAN_A;
    bd.TrigEngSlope1   = TRIGGER_SLOPE_NEGATIVE;
    //bd.TrigEngLevel1   = 128-15;  // negative going pulses
    //bd.TrigEngLevel1   = 102;  // negative going pulses
    //bd.TrigEngLevel1   = 110;  // negative going pulses
    //bd.TrigEngLevel1   = 121;  // negative going pulses: = -5.5 mV for +/- 100 mV full scale
    //bd.TrigEngLevel1   = 118;  // negative going pulses: = -7.8 mV for +/- 100 mV full scale
    bd.TrigEngLevel1   = 116;  // negative going pulses: = -9.4 mV for +/- 100 mV full scale
    //bd.TrigEngLevel1   = 98;  // negative going pulses: = -6 mV for +/- 40 mV full scale
    //bd.TrigEngLevel1   = 128+15;  // positive going pulses
    //bd.TrigEngLevel1   = 128+127;  // positive going pulses
    bd.TriggerEngine2  = TRIG_ENGINE_K;
    bd.TrigEngSource2  = TRIG_DISABLE;
    bd.TrigEngSlope2   = TRIGGER_SLOPE_POSITIVE;
    bd.TrigEngLevel2   = 128;
  } else {
    assert(!"don't know how to generate requested BoardDef");
  }

  return bd;
}

void 
ScopeAlazarATS860::setScopeHandle(HANDLE h) {_h = h;}

HANDLE 
ScopeAlazarATS860::getScopeHandle() { return _h; }

void 
ScopeAlazarATS860::verifyBoardsExist() {

  // ensure that there are boards installed
  U32 nAlazarBoards = AlazarBoardsFound();
  if (nAlazarBoards < 1) {
    cout << "error: no alazar boards found";
    assert(!"error:  no alazar boards found");
  }

  cout << "Found " << nAlazarBoards << " board(s)" << endl;

}

void
ScopeAlazarATS860::verifyBoardType(HANDLE hh) {
  // verify that this board is really an ATS860
  if ( AlazarGetBoardKind(hh) != ATS860 ) {
    cout << "Expecting an ATS860 board (" << ATS860 << ") but found: " 
	 << AlazarGetBoardKind(hh) << endl;
    assert(!"Error: mismatched board type");
  } 
  cout << "Board type verified to be Alazar ATS860 (" << ATS860 
       << ")" << endl;
}

char* ScopeAlazarATS860::GetName() { return "ScopeAlazarATS860"; }

/*
void 
ScopeAlazarATS860::makeScopeConfig(const int sc) {
// Create scope configuration

  _state=sc;

  cout << "ScopeAlazarATS860.makeScopeConfig()" << endl;
  cout << "  sc = " << sc << endl;

  BoardDef bd;

  // Establish some default configuration settings                                
  // to be overridden in the particular cases below.                              
  bd.ClockSource = INTERNAL_CLOCK;

  //switch (sc) {
  // case Scope::SCOPE_SETUP_TEST:
  if (sc == Scope::SCOPE_SETUP_TEST) {
    //case 0:
    // Acquisition info
    bd.RecordCount     = 1;
    //bd.RecLength       = 4*1024;
    //bd.PreDepth        = 0;  // how much pre-trigger data?
    bd.RecLength       = 80*1024;
    bd.PreDepth        = 100*128;  // how much pre-trigger data?
    bd.ClockSource = INTERNAL_CLOCK;      // defaulted
    bd.ClockEdge       = CLOCK_EDGE_RISING;
    //bd.SampleRate      = SAMPLE_RATE_250MSPS;
    bd.SampleRate      = SAMPLE_RATE_25MSPS;

    // Channel A
    bd.CouplingChanA   = DC_COUPLING;
    bd.InputRangeChanA = INPUT_RANGE_PM_10_V;
    //bd.InputRangeChanA = INPUT_RANGE_PM_1_V;
    //bd.InputRangeChanA = INPUT_RANGE_PM_400_MV;
    bd.InputImpedChanA = IMPEDANCE_1M_OHM;
    //bd.InputImpedChanA = IMPEDANCE_50_OHM;
    // Channel B
    bd.CouplingChanB   = DC_COUPLING;
    bd.InputRangeChanB = INPUT_RANGE_PM_1_V;
    bd.InputImpedChanB = IMPEDANCE_50_OHM;
    // Trigger
    //bd.TriEngOperation = TRIG_ENGINE_OP_J_OR_K;
    //bd.TriEngOperation = TRIG_ENGINE_OP_K;
    bd.TriEngOperation = TRIG_ENGINE_OP_J;
    bd.TriggerEngine1  = TRIG_ENGINE_J;
    bd.TrigEngSource1  = TRIG_CHAN_A;
    bd.TrigEngSlope1   = TRIGGER_SLOPE_POSITIVE;
    //bd.TrigEngLevel1   = 160;
    bd.TrigEngLevel1   = 15+128; // ~35% (0 = -100%, 128=0%, 255=100%)
    //bd.TrigEngLevel1   = 128-44; // ~35% (0 = -100%, 128=0%, 255=100%)
    bd.TriggerEngine2  = TRIG_ENGINE_K;
    bd.TrigEngSource2  = TRIG_DISABLE;
    bd.TrigEngSlope2   = TRIGGER_SLOPE_POSITIVE;
    bd.TrigEngLevel2   = 128;
  } else {
    cout << "you should never get here" << endl;
  }

  setScopeConfig(bd);
}
*/

/*
void
ScopeAlazarATS860::clearData() {
  clearCh1Data();
  clearCh2Data();
  cout << "*****************************************" << endl;
  cout << "should also clear the ch1 and ch2 headers" << endl;
  cout << "and the TScopeRecords and ScopeRecords   " << endl;
  cout << "*****************************************" << endl;
}

void 
ScopeAlazarATS860::clearCh1Data() {
  for (vector<ScopeWaveform*>::iterator i=_ch1Data.begin();
       i != _ch1Data.end(); i++) 
    delete *i;
}

void 
ScopeAlazarATS860::clearCh2Data() {
  for (vector<ScopeWaveform*>::iterator i=_ch2Data.begin();
       i != _ch2Data.end(); i++) 
    delete *i;
}
*/
/*
vector<ScopeWaveform*>*
ScopeAlazarATS860::getChannelData(int chNumber) {
  if (chNumber == 1) {
    return &_ch1Data;
  } else if (chNumber == 2) {
    return &_ch2Data;
  } else {
    cout << "getChannelData() Error, unrecognized channel number: " 
	 << chNumber << endl;
    // need a return here...                                                     
    return NULL;                                                               
  }
}
*/
/*
TH1F*
ScopeAlazarATS860::histoOfWaveform(int iwf) {
  // returns a TH1F histogram made from ScopeWaveform data
  TString name = "charge";
  name += iwf;

  BoardDef bd = getScopeConfig();

  TH1F* histo = new TH1F(name, "", bd.RecLength, 0, bd.RecLength);

  for (int ii=0; ii<bd.RecLength; ii++) {
  //for (int ii=0; ii<10; ii++) {
    histo->SetBinContent(ii, (float)(_ch1Data[iwf]->getTrace())[ii]);
  }

  return histo;
}
*/
/*
TScopeRecord*  // made up of ScopeWaveformHeader and TH1F
ScopeAlazarATS860::recordOfWaveform(int iwf) {
  cout << "converting waveform to record: " << iwf << endl;

  TScopeRecord* tsr;
  tsr = new TScopeRecord( histoOfWaveform(iwf), _ch1Hdrs[iwf] );

  return tsr;
}
*/

/*
vector<ScopeRecord*>*
ScopeAlazarATS860::getRecords(int chNumber) {
  if (chNumber == 1) {
    return &_ch1Records;
  } else if (chNumber==2) {
    return &_ch2Records;
  } else {
    cout << "getRecords() Error, unrecognized channel number: "
	 << chNumber << endl;
    return NULL;
  }
}
*/

// the AutoDMA mode allows data to be acquired
// even as previously triggered data is transferred
// from the PCI card to the computer.  This should allow
// for low (no?) dead-time.  The AutoDMA functionality
// is an option for the ATS860, but it is an option that
// we currently do not have...
//int
//ScopeAlazarATS860::acquireTriggersAutoDMA() {
//  cout << "ScopeAlazarATS860::acquireTriggersAutoDMA()" << endl;
//
//  U32 transferLength;
//  BoardDef bd = getScopeConfig();
//  transferLength = bd.RecLength;
//  U32 r1, r2, r3, r4;
//  AUTODMA_STATUS error;
//  HANDLE hh = getScopeHandle();
//  U32 controlFlag = 0;  // start the acquisition right away
//  //U32 controlFlag = 1;  // start the acquisition on next call to AlazarStartCapture
//
//  cout << "AlazarStartAutoDMA()" << endl;
//  status = AlazarStartAutoDMA(hh,
//			      //DMADataBuffers[0],
//			      DataBuffers[0],
//			      UseHeader, channelMode, 
//			      -(long)bd.PreDepth,
//			      transferLength,
//			      RecsPerBuffer,
//			      bd.RecordCount,
//			      &error, controlFlag,
//			      r2, &r3, &r4);
//  if (error != ADMA_Success) {
//    if (error == ADMA_InvalidRecsPerBuffer) {
//      cout << "Error:  RecordCount must be a multiple of RecsPerBuffer\n" 
//	   << "        RecsPerBuffer = " << RecsPerBuffer << endl;
//      return -1;
//    }
//  } else {
//    cout << "*** ADMA_Success ***" << endl;
//  }
//  assert(apiError(status)==0);
//
//  U32 loop_count = 0;
//  bool looping = true;
//  long RecsTransferred = 0;
//  long WhichOne;
//  long TriggersOccurred = 0;
//  int errorCount = 0;
//  int returnValue = 0;
//
//  while (looping) {
//    status = AlazarGetNextAutoDMABuffer(hh, 
//					//DMADataBuffers[0], 
//					//DMADataBuffers[1],
//					DataBuffers[0], 
//					DataBuffers[1],
//					&WhichOne, &RecsTransferred, 
//					&error, r1, r2, 
//					&TriggersOccurred, &r4);
//
//    cout << "status = " << status << endl;
//    cout << "error, ADMA_Success  = " << error << ", " << ADMA_Success << endl;
//
//    // -------------------------------------------
//    // ApiFailed is returned when all of the 
//    // Records per buffer have been transferred
//    if (status == ApiFailed) {
//      loop_count = 0;
//      looping = false;
//    }
//
//    // ----------------------------------
//    // Valid Data exists in either of the buffers
//    if ( (WhichOne == 0) || (WhichOne == 1) ) {
//      cout << "VALID DATA " << endl;
//      //SaveToChannelFiles(DMADataBuffers[WhichOne]);
//      //processBufferData(DMADataBuffers[WhichOne]);
//      processBufferData(DataBuffers[WhichOne]);
//      loop_count = 0;
//    }
//
//    // ---------------------------------
//    // an in-elegant software delay
//    loop_count += 1;
//    //if (loop_count > 1000000) {
//    if (loop_count > 10) {
//      cout << " ******** loop_count exceeded max value *********" << endl;
//      cout << "TriggersOccurred = " << TriggersOccurred << endl;
//      if ( (TriggersOccurred > 0) && (TriggersOccurred < (long)RecsPerBuffer) ) {
//	//AlazarAbortAutoDMA(hh, DMADataBuffers[0], &error, r1, r2, &r3, &r4);
//	AlazarAbortAutoDMA(hh, DataBuffers[0], &error, r1, r2, &r3, &r4);
//      }
//      looping = false;
//      cout << " ******** software timeout exceeded *********" << endl;
//      errorCount += 1;
//    }
//
//    // ---------------------------------
//    // An AutoDMA error has happened, so abort the current acquisition
//    if (error == ADMA_OverFlow) {
//      looping = false;
//      returnValue = -4;
//    }
//  }
//  AlazarCloseAUTODma(hh);
//  
//  // close the file handles
//  
//  // free the memory
//  freeMemory();
//
//  return returnValue;
//}

/*
int 
ScopeAlazarATS860::processBufferData(U8* Data) {
  //DMADataBuffers[WhichOne]);

  //channelMode = CHANNEL_A;
  if (channelMode == (CHANNEL_A | CHANNEL_B)) {
    cout << "processBufferData(): not yet implemented" << endl;
  } else {
    cout << "processBufferData() error: unknown channelMode = " << channelMode << endl
	 << "CHANNEL_A = " << CHANNEL_A << ", CHANNEL_B = " << CHANNEL_B << endl;
    return -1;
  }

  return 0;
}
*/


void 
ScopeAlazarATS860::freeMemory() {
  //cout << "need to do memory clean-up here" << endl;
}

int 
ScopeAlazarATS860::apiError(RETURN_CODE status) {
  if (status != ApiSuccess) {
    freeMemory();
    cout << "************************************" << endl;
    cout << "  Alazar API failure -- status = " << status << endl;
    cout << "************************************" << endl;
    return 1;
  } else {
    //cout << "   alazar API success" << endl;
  }
  return 0;
}

void 
ScopeAlazarATS860::printDriverVersion() {
  U8 drvVerMaj, drvVerMin, drvVerRev; // major, minor and revision numbers
  AlazarGetDriverVersion(&drvVerMaj, &drvVerMin, &drvVerRev);
  cout << "Alazar driver version = " << (int)drvVerMaj << "." 
       << (int)drvVerMin << "." << (int)drvVerRev << endl;
}

void 
ScopeAlazarATS860::printSDKVersion() {
  U8 sdkVerMaj, sdkVerMin, sdkVerRev;
  AlazarGetSDKVersion(&sdkVerMaj, &sdkVerMin, &sdkVerRev);
  cout << "Alazar SDK version    = " << (int)sdkVerMaj << "." 
       << (int)sdkVerMin << "." << (int)sdkVerRev << endl;
}

void
ScopeAlazarATS860::printScopeMemorySize(HANDLE hh) {
  U32 retVal, dummy=0;
  AlazarQueryCapability( hh, MEMORY_SIZE, dummy, &retVal );
  cout << "  This board has a memory size of " << (int)retVal 
       << " samples per channel" << endl;
}

void
ScopeAlazarATS860::printScopeChannelInfo(HANDLE hh) {
  U32 MemSize;
  U8 SampleSize;
  AlazarGetChannelInfo(hh, &MemSize, &SampleSize);
  cout << "Acquisition memory size per channel: " << (int)MemSize << endl;
  cout << "Number of bits per sample: " << (int)SampleSize << endl;
}


void
ScopeAlazarATS860::printScopeSerialNumber(HANDLE hh) {
  // GET_SERIAL_NUMBER is not defined in SDK 4.9.1 (it is in for 5.X)
  //U32 retVal, dummy=0;
  //AlazarQueryCapability( hh, GET_SERIAL_NUMBER, dummy, &retVal );
  //cout << "  Serial number = " << (int)retVal << endl;
}

void 
ScopeAlazarATS860::printScopeSetup() {
  ScopeConfig* sc = getScopeConfig();

  cout << GetName() << ":  Scope setup" << endl;
  cout << "  Sampling Rate = " << getSamplingRate()*1e-6 << " MS/s" << endl;
  cout << "  timespan      = " << 1e6*sc->board()->getRecordSize()/getSamplingRate() << " us" << endl;
  U32 level = sc->trig()->trigLevel1();
  cout << "  trigger level = " << voltsOfLevel(0, level) << " volts (" << level << " unitless)" << endl;

  for (int ii=0; ii<getNChannels(); ii++) {
    cout << "  Channel " << ii << endl;
    cout << "    Voltage Range = " << getVoltageMin(ii) << ", " 
	 << getVoltageMax(ii) << endl;
  }
}


/*
  // Allocate the data buffers for AutoDMA data transfer.
  } else if (_acqMode == ACQ_MODE_AUTODMA) {
    cout << "_acqMode == ACQ_MODE_AUTODMA" << endl;
    bytesPerBuffer = RecsPerBuffer * (1+(channelMode==(CHANNEL_A | CHANNEL_B)))
      * (bd.RecLength * sizeof(U8) + sizeof(ALAZAR_HEADER)*(UseHeader != 0));
    cout << "bytesPerBuffer = " << bytesPerBuffer << endl;
*/

