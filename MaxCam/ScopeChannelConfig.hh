// ScopeChannelConfig.hh
//
#ifndef SCOPE_CHANNEL_CONFIG_HH
#define SCOPE_CHANNEL_CONFIG_HH

//#include <string>

#include "TObject.h"
#include "TString.h"

#include "ScopeTypes.hh"

//#include "AlazarApi.h"
//#include "AlazarCmd.h"

//
// Class to contain the configuration for the scope channel
// 
// 
class ScopeChannelConfig : public TObject {

public:

  // Destructor
  ~ScopeChannelConfig();
  // Constructors
  ScopeChannelConfig();
  ScopeChannelConfig(AU8  chanNum, AU32 coupling, 
		     AU32 vRange,  AU32 impedance, TString* label);
		     //AU32 vRange,  AU32 impedance, std::string* label);

  //std::string* label();
  TString* label();

  AU8  getChannelNumber();
  AU32 getCoupling();
  AU32 getVoltageRange();
  AU32 getInputImpedance();

  void setVoltageRange(AU32 range);
  void setInputCoupling(AU32 coupling);
  void setInputImpedance(AU32 impedance);

private:
  AU8  _channelNumber;
  AU32 _coupling;
  AU32 _voltageRange;
  AU32 _impedance;
  //std::string* _label;
  TString* _label;
  char* GetName();

ClassDef(ScopeChannelConfig, 1)
};

#endif
