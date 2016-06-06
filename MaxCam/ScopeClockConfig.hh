// ScopeClockConfig.hh
//
#ifndef SCOPE_CLOCK_CONFIG_HH
#define SCOPE_CLOCK_CONFIG_HH

#include "TObject.h"

//#include "AlazarApi.h"
//#include "AlazarCmd.h"

#include "ScopeTypes.hh"

//
// Class to contain the configuration for the scope clock
// 
// 
class ScopeClockConfig : public TObject {

public:

  // Destructor
  ~ScopeClockConfig();
  // Constructors
  ScopeClockConfig();
  ScopeClockConfig(AU32 cs, AU32 cr, AU32 ce, AU32 cd);

  void setClockSource(AU32 cs);
  void setClockRate(AU32 cr);
  void setClockEdge(AU32 ce);
  void setClockDecimation(AU32 cd);

  AU32 source();

  AU32 getClockSource();
  AU32 getClockRate();
  AU32 getClockEdge();
  AU32 getClockDecimation();

private:
  AU32 _clockSource;
  AU32 _clockRate;
  AU32 _clockEdge;
  AU32 _clockDecimation;

ClassDef(ScopeClockConfig, 1)

};

#endif
