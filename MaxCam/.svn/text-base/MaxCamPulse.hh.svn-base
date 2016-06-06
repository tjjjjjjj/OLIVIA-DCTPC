#ifndef MAXCAM_PULSE_HH
#define MAXCAM_PULSE_HH

#include "TH1F.h"
#include <ostream>
#include <iostream>

class MaxCamPulse : public TObject
{

  public: 
    MaxCamPulse(){init = false;}
    ~MaxCamPulse();
    MaxCamPulse(Int_t nbin);

    /** the bin in which the pulse maximum occurs */
    Int_t getBin() { return _nbin; }
  
    /** calculated as the maximum of pulse, minus the baseline*/
    void setPulseHeight(Double_t value) { _pulseHeight=value; }
    Double_t getPulseHeight() { return _pulseHeight; }

    /** calculated as the center of the time bin in which the pulse height occurs*/
    void setPulseHeightTime(Double_t value) { _pulseHeightTime=value; }
    Double_t getPulseHeightTime() { return _pulseHeightTime; }

    // time center of first bin in pulse 
    void setPulseStartTime(Double_t value) { _pulseStartTime=value; }
    Double_t getPulseStartTime() { return _pulseStartTime; }

    // first bin of pulse 
    void setPulseStartBin(Double_t value) { _pulseStartBin=value; }
    Double_t getPulseStartBin() { return _pulseStartBin; }

    // time center of last bin of pulse 
    void setPulseEndTime(Double_t value) { _pulseEndTime=value; }
    Double_t getPulseEndTime() { return _pulseEndTime; }

    // last bin of pulse 
    void setPulseEndBin(Double_t value) { _pulseEndBin=value; }
    Double_t getPulseEndBin() { return _pulseEndBin; }

    // integral of pulse 
    void setPulseIntegral(Double_t value) { _pulseIntegral=value; }
    Double_t getPulseIntegral() { return _pulseIntegral; }

    /** prints quantities for the pulse */
    void print(ostream & out = std::cout);

    private:

      bool init; 

      Int_t _nbin;
      Double_t _pulseHeight;
      Double_t _pulseHeightTime;
      Double_t _pulseStartTime;
      Double_t _pulseStartBin;
      Double_t _pulseEndTime;
      Double_t _pulseEndBin;
      Double_t _pulseIntegral;
   
    ClassDef(MaxCamPulse,1); 
    
};

#endif

