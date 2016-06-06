#include "MaxCamTriggerGroup.hh"
#include <strings.h>

ClassImp(MaxCamTriggerGroup);


MaxCamTriggerGroup::~MaxCamTriggerGroup()
{
  _waveforms->SetOwner(kTRUE); 
  _waveforms->Clear(); 
  delete _waveforms; 
}


MaxCamTriggerGroup::MaxCamTriggerGroup()
{
  _waveforms = new TObjArray(); 
}

MaxCamTriggerGroup::MaxCamTriggerGroup(unsigned int i)
{
  _waveforms = new TObjArray(); 
  _idx = i; 
}

MaxCamWaveform * MaxCamTriggerGroup::getWaveform(char * type)
{

  for (int i = 0; i < nWaveForms(); i++)
  {
    if (strcmp(getWaveform(i)->getType(),type)==0)
    {
      return getWaveform(i);
    }
  }

  return NULL; 
}

