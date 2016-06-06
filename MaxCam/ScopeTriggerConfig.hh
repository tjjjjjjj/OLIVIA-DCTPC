// ScopeTriggerConfig.hh
//
#ifndef SCOPE_TRIGGER_CONFIG_HH
#define SCOPE_TRIGGER_CONFIG_HH

#include "TObject.h"

//#include "AlazarApi.h"
//#include "AlazarCmd.h"
#include "ScopeTypes.hh"

//
// Class to contain the configuration for the scope trigger
// 
// should probably make a new class called ScopeTriggerEngineConfig()
// which serve as data members of this class.  There would be
// one ScopeTriggerEngineConfig() per trigger engine (the alazar board
// has two trigger engines).  Any general trigger information, such 
// as TriggerOperation -- see _trgOp below -- 
// (e.g. J, K, J_OR_K, J_AND_OR_K, etc) would live in this class alone.
// 
class ScopeTriggerConfig : public TObject {

public:

  // Destructor
  ~ScopeTriggerConfig();
  // Constructors
  ScopeTriggerConfig();
  ScopeTriggerConfig(AU32 trgOp, 
		     AU32 trgEn1, AU32 trgSrc1, AU32 trgSlope1, AU32 trgLevel1,
		     AU32 trgEn2, AU32 trgSrc2, AU32 trgSlope2, AU32 trgLevel2);

  // get
  AU32 trigOp();
  AU32 trigEn1();
  AU32 trigEn2();
  AU32 trigSrc1();
  AU32 trigSrc2();
  AU32 trigSlope1();
  AU32 trigSlope2();
  AU32 trigLevel1();
  AU32 trigLevel2();

  // set
  void setTrigOp(AU32 to);

  void setTrigLevel1(AU32 tl);
  void setTrigEngine1(AU32 te);
  void setTrigSource1(AU32 ts);
  void setTrigSlope1(AU32 ts);

  void setTrigLevel2(AU32 tl);
  void setTrigEngine2(AU32 te);
  void setTrigSource2(AU32 ts);
  void setTrigSlope2(AU32 ts);
  
private:

  AU32 _trgOp;
  AU32 _trgEn1;
  AU32 _trgSrc1;
  AU32 _trgSlope1;
  AU32 _trgLevel1;
  AU32 _trgEn2;
  AU32 _trgSrc2;
  AU32 _trgSlope2;
  AU32 _trgLevel2;

ClassDef(ScopeTriggerConfig, 1)

};

#endif
