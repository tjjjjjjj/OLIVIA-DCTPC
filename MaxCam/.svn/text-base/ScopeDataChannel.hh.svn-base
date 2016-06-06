//ScopeDataChannel.hh
//
#ifndef SCOPE_DATA_CHANNEL_HH
#define SCOPE_DATA_CHANNEL_HH

#include <vector>

#include "TObject.h"
#include "TString.h"

class ScopeWaveformData;

//
// Class to hold a data for a single channel of a scope.a
// Has data member vector<ScopeWaveformData*> that holds the data
// read out from that channel on a single trigger.
//
class ScopeDataChannel : public TObject {

public:

  // Destructor
  virtual ~ScopeDataChannel();
  // Constructors
  ScopeDataChannel();
  //ScopeDataChannel(int chNum, unsigned int recLength, 
  //int bytesPerBuffer, std::string* chLabel);
  //ScopeDataChannel(int chNum, std::string* chLabel, int bytesPerBuffer);
  ScopeDataChannel(int chNum, TString* chLabel, int bytesPerBuffer);
    
  TString* label();
  int chNum();

  int addWaveform( ScopeWaveformData* swf );
  ScopeWaveformData* wf(int id);

  // should really be U8
  unsigned char* buffer();

  //int getBytesPerBuffer();
  //unsigned int getRecLength();

  std::vector<ScopeWaveformData*>* wfs();
  std::vector<float>* wftimes();

  void addTime(float seconds) { _wfTime.push_back(seconds); }
  float getTime(int i) { return _wfTime[i]; }
    
private:
  std::vector<ScopeWaveformData*> _waveforms;
  int    _channelNumber;
  //should really be U32
  //unsigned int    _recLength;
  //int    _bytesPerBuffer;
  TString _channelLabel;

  // A pointer to the first position of an array which is 
  // somewhere in memory.  The scope will use this pointer
  // to put data into the array when it transfers data from
  // the board to memory.
  // Should really be U8*
  unsigned char* _dataBufferPtr;  

  void allocateMemoryForBuffer(int bytesPerBuffer);

  char* GetName();

  std::vector<float> _wfTime;
    
  ClassDef(ScopeDataChannel,0)
};

#endif
