//ScopeDataChannel.cc
//

//#include <string>
#include <stdlib.h>

#include <iostream>
#include <vector>

#include "AlazarApi.h"
#include "AlazarCmd.h"

#include "ScopeDataChannel.hh"
//#include "ScopeConfig.hh"
#include "ScopeWaveformData.hh"
#include "TROOT.h"


using std::cout;
using std::endl;
//using std::string;
using std::vector;

ClassImp(ScopeDataChannel)

ScopeDataChannel::~ScopeDataChannel() {
  cout << GetName() << ": destructor" << endl;
  
  // remove all child waveforms
  for (vector<ScopeWaveformData*>::iterator ii=_waveforms.begin();
       ii != _waveforms.end(); ++ii) {
    gROOT->Delete((*ii)->GetName());
    delete *ii;
  }
  
  _wfTime.clear();

  //delete _dataBufferPtr;
}

ScopeDataChannel::ScopeDataChannel() {}

//ScopeDataChannel::ScopeDataChannel(int chNum, string* chLabel, int bytesPerBuffer) {
ScopeDataChannel::ScopeDataChannel(int chNum, TString* chLabel, int bytesPerBuffer) {
  //The ATS860 is an 8-bit device.
  //In traditional mode, you must pad by 16 samples.
  //You therefore need:
  //   bytesPerBuffer = (bd.RecLength+16)*sizeof(U8)
  cout << GetName() << ": constructor" << endl;
  _channelNumber  = chNum;
  _channelLabel   = *chLabel;
  _dataBufferPtr  = NULL;
  allocateMemoryForBuffer(bytesPerBuffer);
  //_recLength      = recLength;
  //_bytesPerBuffer = bytesPerBuffer;
  //allocateMemoryForBuffer(_bytesPerBuffer);
}

//string* 
TString* 
ScopeDataChannel::label() { return &_channelLabel; }

int
ScopeDataChannel::chNum() { return _channelNumber; }

int
ScopeDataChannel::addWaveform( ScopeWaveformData* swf ) {
  _waveforms.push_back(swf); 
  return _waveforms.size() - 1;
}
ScopeWaveformData* 
ScopeDataChannel::wf(int ii) {
  // get the ii-th waveform from the _waveforms vector
  return _waveforms.at(ii);
}

vector<ScopeWaveformData*>*
ScopeDataChannel::wfs() { return &_waveforms; }

vector<float>*
ScopeDataChannel::wftimes() { return &_wfTime; }

U8* 
ScopeDataChannel::buffer() { return _dataBufferPtr; }

//int 
//ScopeDataChannel::getBytesPerBuffer() { return _bytesPerBuffer; }

//U32 
//ScopeDataChannel::getRecLength() { return _recLength; }

//
// PRIVATE METHODS
//
char* 
ScopeDataChannel::GetName() { return "ScopeDataChannel"; }

void 
ScopeDataChannel::allocateMemoryForBuffer(int bytesPerBuffer) {
  _dataBufferPtr = (U8 *) malloc(bytesPerBuffer);
}

