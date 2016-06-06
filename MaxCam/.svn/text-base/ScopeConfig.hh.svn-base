//ScopeConfig.hh
//
#ifndef SCOPE_CONFIG_HH
#define SCOPE_CONFIG_HH

#include <vector>

#include "TObject.h"

//#include "AlazarApi.h"
//#include "AlazarCmd.h"

#include "ScopeTypes.hh"

#include "ScopeBoardInfo.hh"
#include "ScopeClockConfig.hh"
#include "ScopeTriggerConfig.hh"
#include "ScopeChannelConfig.hh"
//
// class to hold scope configuration information
//
class ScopeConfig : public TObject {

public:

  // Destructor
  ~ScopeConfig();
  // Constructors
  ScopeConfig();

  // set up records, clock and trigger, but not channels 
  // (use addChannelConfig() to add channels)
  ScopeConfig(ScopeBoardInfo* boardcfg,
	      ScopeClockConfig* clockcfg,
	      ScopeTriggerConfig* trgcfg);
  
  // set up a 2 channel scope
  ScopeConfig(ScopeBoardInfo* boardcfg,
	      ScopeClockConfig* clockcfg,
	      ScopeTriggerConfig* trgcfg, 
	      ScopeChannelConfig* chancfg1, 
	      ScopeChannelConfig* chancfg2);

  // copy constructor
  ScopeConfig(ScopeConfig& sc);

  ScopeConfig(int junk);

  void setBoardConfig(ScopeBoardInfo* bi);
  void setClockConfig(ScopeClockConfig* cfg);
  void setTriggerConfig(ScopeTriggerConfig* cfg);
  void addChannelConfig(ScopeChannelConfig* cfg);

  ScopeBoardInfo*     board();
  ScopeClockConfig*   clock();
  ScopeTriggerConfig* trig();
  ScopeChannelConfig* chan(int ii);
  std::vector<ScopeChannelConfig*> chans();

  void print();
  void printall();

    virtual const char* GetName() const;

private:
  ScopeBoardInfo*     _boardCfg;
  ScopeClockConfig*   _clockCfg;
  ScopeTriggerConfig* _trgCfg;
  std::vector<ScopeChannelConfig*> _chanCfg;

ClassDef(ScopeConfig, 1)
};

#endif
