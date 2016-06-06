#include "MaxCamWaveform.hh"
#include <strings.h>

ClassImp(MaxCamWaveform)

MaxCamWaveform::~MaxCamWaveform()
{
  if (init)
  {
   _wf->Delete(); 
  }
}


MaxCamWaveform::MaxCamWaveform(TH1F * h, unsigned int trigger_group, const char * type, bool invert, int nbaselinebins, int nrebin)
{
  _wf = (TH1F*) h->Clone();  
  MaxCamWaveformTools * tools = new MaxCamWaveformTools(_wf,invert,nbaselinebins,nrebin); 

  _baseline = tools->getBaseline();
  _baselineRMS = tools->getBaselineRMS();
  _peakHeight = tools->getPeakHeight();
  _peakHeightTime = tools->getPeakHeightTime();
  _riseTime = tools->getRiseTime();
  _fallTime = tools->getFallTime();;
  _troughDepth = tools->getTroughDepth();
  _wfMinimum = tools->getWfMinimum();
  _peakAveragePosTime = tools->getPeakAveragePosTime();
  _peakAverageNegTime = tools->getPeakAverageNegTime();

  _trigger_group = trigger_group; 
  _type = TString(type); 
  init = true; 

  delete tools; 
}


void MaxCamWaveform::print(ostream & out)
{
  out << "MaxCamWaveform of type " << _type << " in trigger group " << _trigger_group << std::endl; 
  out << "baseline: " << getBaseline() << " baselineRMS: " << getBaselineRMS() << " peakHeight: " << getPeakHeight() << std::endl; 
  out << "peakHeightTime: " << getPeakHeightTime() << " riseTime: " << getRiseTime() << " fallTime: " << getFallTime() << std::endl; 
  out << "troughDepth: " << getTroughDepth() << " wfMinimum: " << getWfMinimum() << std::endl; 
  out  << "peakAveragePosTime: " << getPeakAveragePosTime() << " peakAverageNegTime: " << getPeakAverageNegTime() << std::endl; 
}
