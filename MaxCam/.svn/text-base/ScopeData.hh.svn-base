// class to hold a scope data
//

#ifndef SCOPE_DATA_HH
#define SCOPE_DATA_HH

#include <vector>

#include "ScopeConfig.hh"
#include "ScopeDataChannel.hh"
#include "ScopeWaveformData.hh"
//class ScopeDataChannel;
//class ScopeWaveformData;

class ScopeData {

public:
  // Destructor
  ~ScopeData();
  // Constructors
  ScopeData();
  ScopeData(ScopeConfig* sc);

  ScopeConfig* config();
  void setScopeConfig(ScopeConfig* sc);
  ScopeConfig* getScopeConfig();

  ScopeDataChannel* dataChan(int id);
  int addDataChannel(ScopeDataChannel* dataChan);

  std::vector<ScopeWaveformData*> getWaveforms();

  int  getNValidTriggers();
  void setNValidTriggers(int nval);

  void clearData();

private:
  int _nValidTriggers;

  std::vector<ScopeDataChannel*> _channels;
  ScopeConfig* _scopeConfig;

  void initDataset();
  char* GetName();

};

#endif
