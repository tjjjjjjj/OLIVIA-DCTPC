#ifndef MAXCAM_TRIGGER_GROUP_HH
#define MAXCAM_TRIGGER_GROUP_HH
#include "MaxCamWaveform.hh"
#include "TObjArray.h"
#include <vector>


class MaxCamTriggerGroup : public TObject
{

  public: 

    MaxCamTriggerGroup();
    MaxCamTriggerGroup(unsigned int idx);
    ~MaxCamTriggerGroup(); 
    
    TObjArray * getWaveForms(){return _waveforms;}
    void addWaveform(MaxCamWaveform* wf){_waveforms->Add(wf);}
    unsigned int nWaveForms() { return _waveforms->GetEntries(); } 
    MaxCamWaveform * getWaveform(unsigned int index){return (MaxCamWaveform*) _waveforms->At(index);}     
    MaxCamWaveform * getWaveform(char * type); 
    unsigned int getIndex(){ return _idx;} 

  private: 
     
    TObjArray * _waveforms; 
    unsigned int _idx; 
 
  ClassDef(MaxCamTriggerGroup,1);
};

#endif 
