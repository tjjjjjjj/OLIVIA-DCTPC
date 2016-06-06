#include "cleanSkimProcessor.hh"
#include "../../../MaxCam/DmtpcDataset.hh"
#include "../../../MaxCam/DmtpcEvent.hh"
#include "../../../MaxCam/ScopeDataInfo.hh"

CleanSkimConfig * skim::preprocess(DmtpcDataset * d, const char * seed)
{

  gROOT->cd();
  CleanSkimConfig *cfg = new CleanSkimConfig(seed); 
  gROOT->cd();
  d->getEvent(0); 


  //Check for overscan  
  if (d->event()->overscan() == 0 || d->getBiasFrameOverscan(1) == 0)  
  {
    cfg->setHasOverscan(false); 
  }
  
  //Check for charge readout 
  if (d->event()->scopeData() == 0) 
  {
    cfg->setNScopeChannels(0); 
  }

  // only works if there's a trigger in the first event
  /*
  else if (d->event()->scopeDataInfo()) 
  {
    //If we have charge readout, get the number of channels
    int nch = d->event()->scopeData()->GetEntries(); 
    if (nch > 0) 
      nch /= ((ScopeDataInfo*)d->event()->scopeDataInfo()->At(0))->getNTriggers(); 
    cfg->setNScopeChannels(nch); 
  }
  */

  return cfg; 
}







