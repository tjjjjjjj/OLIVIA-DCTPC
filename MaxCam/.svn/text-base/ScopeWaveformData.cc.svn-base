// ScopeWaveformData.cc
//
#include "waveformtools/include/WaveformTools.hh" 
#include "waveformtools/include/WaveformAnalysis.hh" 
#include "DmtpcRootTools.hh"
#include <iostream>

#include "ScopeWaveformData.hh"

using std::cout;
using std::endl;

ClassImp(ScopeWaveformData)

ScopeWaveformData::~ScopeWaveformData() {
  //cout << GetName() << ": destructor" << endl;
}
ScopeWaveformData::ScopeWaveformData() { 
  _signed = false;
  _nresample = 1;
}
//ScopeWaveformData::ScopeWaveformData(const ScopeWaveformData &other) {
//  _nSamples    = other._nSamples;
//  _dataChannel = other.sdc;
//  cout << GetName() << ": copy constructor not done" << endl;
//}

ScopeWaveformData::ScopeWaveformData(TString histoName, int nSamples, double timeStamp, float leveltoVolt, short zero) : TH1D(histoName, "", nSamples, 0, nSamples) {
  _nSamples = nSamples;
  _timeStamp = timeStamp;
  _leveltoVolt = leveltoVolt;
  _zero = zero;
  _signed = false; 
  _nresample = 1; 

}

ScopeWaveformData::ScopeWaveformData(TString histoName, int nSamples, float xlow, float xup, double timeStamp, float leveltoVolt, short zero) : TH1D(histoName, "", nSamples, xlow, xup) {
  _nSamples = nSamples;
  _nresample = 1; 
  _timeStamp = timeStamp;
  _leveltoVolt = leveltoVolt;  
  _zero = zero;
  _signed = false; 
}

TH1S* ScopeWaveformData::getAsShortHist() const 
{
  TH1S * out = (TH1S*) DmtpcRootTools::newTH1StealSize(this,'S',fName,fTitle); 

  for (int i = 1; i <= GetNbinsX(); i++)
  {
    out->SetBinContent(i, ((unsigned char)GetBinContent(i)) ); 
  }

 return out; 
}

TH1D* ScopeWaveformData::getAsSignedHist() const 
{
  TH1D * out = (TH1D*) DmtpcRootTools::newTH1StealSize(this,'C',fName,fTitle); 

  for (int i = 1; i <= GetNbinsX(); i++)
  {
    out->SetBinContent(i, char(short((unsigned char)GetBinContent(i))-_zero) ); 
  }

 return out; 
}

void ScopeWaveformData::convertToSignedHist()
{
  for (int i = 1; i <= GetNbinsX(); i++)
  {
    SetBinContent(i, char(short((unsigned char)GetBinContent(i))-_zero) ); 
  }
  _signed = true; 
}


int ScopeWaveformData::zeroSuppress(double nsigma, unsigned nrequired, unsigned window_size, int nsamples, double max_deviation )
{

  double noise; 
  double base = waveform::analysis::baseline(this, noise,1,nsamples); 

  std::cout << "base: "  << base << std::endl; 
  std::cout << "noise: " << noise << std::endl; 
  double zero = _signed ? 0 : _zero;
  //baseline isn't close to 0... we shouldn't zero suppress
  if (fabs(base-zero) > max_deviation * noise)
  {
    return 0; 
  }

  _nsuppressed =  waveform::tools::zeroSuppress(this, zero, nsigma * noise, nrequired, window_size,0,&noise); 
  _suppressed_rms = float(noise); 
  return _nsuppressed; 
}


void ScopeWaveformData::resample(double factor)
{
  _nresample = factor; 
  TH1 * tmp = waveform::tools::resample(this,factor); 
  this->SetBins(tmp->GetNbinsX(), tmp->GetXaxis()->GetXmin(), tmp->GetXaxis()->GetXmax()); 
  for (int i = 1; i <= tmp->GetNbinsX(); i++)
  {
    SetBinContent(i,tmp->GetBinContent(i)); 
  }

  delete tmp; 
}

void ScopeWaveformData::cropZeros()
{
  TH1 * tmp = waveform::tools::cropZeros((TH1D*)this); 
  this->SetBins(tmp->GetNbinsX(), tmp->GetXaxis()->GetXmin(), tmp->GetYaxis()->GetXmax()); 
  for (int i = 1; i <= tmp->GetNbinsX(); i++)
  {
    SetBinContent(i,tmp->GetBinContent(i)); 
  }

  delete tmp; 
}
