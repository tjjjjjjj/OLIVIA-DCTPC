// ScopeHandler.hh
//
#ifndef SCOPE_HANDLER_CAEN_HH
#define SCOPE_HANDLER_CAEN_HH

//#include <string>
#include <vector>

#include "TString.h"
#ifdef SCOPE_CAEN
#include <CAENComm.h>
#include <CAENVMEoslib.h>
#include <CAENDigitizer.h>
#else
#include "AlazarApi.h"
#include "AlazarCmd.h"
#endif

#include "ScopeTypes.hh"

class ScopeData;
class ScopeDataInfo;
class ScopeWaveformData;
class ScopeConfig;

//
// Class to run an oscilloscope
//
class ScopeHandlerCAEN {

public:

  // Destructor
  ~ScopeHandlerCAEN();
  // Constructor
  ScopeHandlerCAEN();

private:
  char* GetName();
};


#endif
