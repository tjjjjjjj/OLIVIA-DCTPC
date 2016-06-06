// ScopeWaveform.hh
//
#ifndef SCOPE_WAVEFORM_HH
#define SCOPE_WAVEFORM_HH

#include "TH1F.h"
//#include "TString.h"

#include "ScopeDataChannel.hh"

class TString;

//
// Class to contain waveforms from an oscilloscope
// 
// 
class ScopeWaveform : public TH1F { 

public:

  // Destructor
  ~ScopeWaveform();

  // Constructors
  ScopeWaveform();
  //ScopeWaveform(const ScopeWaveform &other);  // copy constructor
  ScopeWaveform(ScopeDataChannel* sdc, TString histoName, int nSamples);
  
  ScopeWaveform(ScopeDataChannel* sdc, TString histoName, 
		int nSamples, float xlow, float xup);

  ScopeDataChannel* chan();

private:  
  int _nSamples;
  ScopeDataChannel* _dataChannel;
  //char* GetName();
  ClassDef(ScopeWaveform,1)
};

#endif
