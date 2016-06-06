// ScopeWaveform.cc
//
#include <iostream>

#include "ScopeWaveform.hh"

using std::cout;
using std::endl;

ClassImp(ScopeWaveform)

ScopeWaveform::~ScopeWaveform() {
  //cout << GetName() << ": destructor" << endl;
}
ScopeWaveform::ScopeWaveform() {}
//ScopeWaveform::ScopeWaveform(const ScopeWaveform &other) {
//  _nSamples    = other._nSamples;
//  _dataChannel = other.sdc;
//  cout << GetName() << ": copy constructor not done" << endl;
//}

ScopeWaveform::ScopeWaveform(ScopeDataChannel* sdc, TString histoName, int nSamples) : TH1F(histoName, "", nSamples, 0, nSamples) {
  //cout << GetName() << ": constructor" << endl;
  _dataChannel = sdc;
  _nSamples = nSamples;
}

ScopeWaveform::ScopeWaveform(ScopeDataChannel* sdc, TString histoName, int nSamples, float xlow, float xup) : TH1F(histoName, "", nSamples, xlow, xup) {
  //cout << GetName() << ": constructor" << endl;
  _dataChannel = sdc;
  _nSamples = nSamples;
}

ScopeDataChannel* 
ScopeWaveform::chan() { return _dataChannel; }

//
// PRIVATE methods
//
//char* ScopeWaveform::GetName() { return "ScopeWaveform"; }
