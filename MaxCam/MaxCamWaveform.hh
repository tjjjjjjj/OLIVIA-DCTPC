#ifndef MAXCAM_WAVEFORM_HH
#define MAXCAM_WAVEFORM_HH

#include "TH1F.h"
#include "MaxCamWaveformTools.hh" 
#include <ostream>
#include <iostream>

/**
 * This class is intended to store a waveform in a skim tree. It largely mirrors
 * MaxCamWaveformTools
 *
 */

class MaxCamWaveform : public TObject
{

  public: 
    MaxCamWaveform(){init = false;}
    ~MaxCamWaveform();
    MaxCamWaveform(TH1F * wf, unsigned int trigger_group, const char * type, bool invert = false, int nbaselinebins = 25, int nrebin = 1);


    const char * getType(){ return _type;} 
    unsigned int getTriggerGroup() { return _trigger_group; } 
    TH1F * getWaveform() { return _wf; } 

   /** calculated as the average of the first nbaselinebins */
    float getBaseline() { return _baseline; }

   /** calculated as the rms of the first nbaselinebins */
    float getBaselineRMS() { return _baselineRMS; }

   /** calculated as the maximum of TH1, minus the baseline*/
    float getPeakHeight() { return _peakHeight; }

   /** calculated as the center of the time bin in which the peak height occurs*/
    float getPeakHeightTime() { return _peakHeightTime; }

   /** calculated as the minimum of the TH1, minus the baseline */
    float getMinimum() { return _wfMinimum; }

   /** calculated as the time to go from 10% to 90% of the peak height */
    float getRiseTime() { return _riseTime; }

   /** same as minimum */
    float getTroughDepth() { return _troughDepth; }

   /** calculated as the time to fall from 90% to 10% of peak. If the waveform does not get down to 10% of peak, the fall time is calculated as 90% of peak to the end of readout */
    float getFallTime() { return _fallTime; }

   /** same as Minimum */
    float getWfMinimum() { return _wfMinimum; }

   
    float getPeakAveragePosTime() { return _peakAveragePosTime; }
    float getPeakAverageNegTime() { return _peakAverageNegTime; }
    
   /** prints quantities for the waveform */
    void print(ostream & out = std::cout);

    private:
      unsigned int _trigger_group; 
      TH1F * _wf; 
      TString _type; 
      bool init; 
   
      float _baseline;
      float _baselineRMS;
      float _peakHeight;
      float _peakHeightTime;
      float _riseTime;
      float _fallTime;
      float _troughDepth;
      float _wfMinimum;
      float _peakAveragePosTime;
      float _peakAverageNegTime;
   

    ClassDef(MaxCamWaveform,1); 
    
};

#endif

