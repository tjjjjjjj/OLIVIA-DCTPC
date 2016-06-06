//
// LabJackU3.hh
//

#ifndef LABJACK_U3_HH
#define LABJACK_U3_HH

#include "TString.h"

#ifndef __CINT__
#include "u3.h"
#else
struct u3CalibrationInfo;
#endif

#include "labjackusb.h"


//
// Class to control the LabJack U3-LV multifunction DAQ device
// Has analog in, analog out, digital in, timers, etc.
// For now, I just want to control the DAC (analog output).
// 

class LabJackU3 {

public:
  // Destructor
  virtual ~LabJackU3(){;}
  // Constructors
  LabJackU3();
  float setDAC(Float_t volts, Int_t channel=0);
  long getUpTime();

  /** Sets up two counters at pinOffset and pinOffset +1. No timers are enabled.
   *
   *  If in the future, we want to enable timers, we will need something like setUpTimersAndCounters. 
   *  \param pinOffset the pin offset to use for the counters 
   *  \returns 0 on success
   * **/ 
  long setUpCounters(long pinOffset = 8);  

  /** Reads the values of the two counters. Pass pointers to the values you want to update or 0 if you don't want to read 
   *  that counter.
   *  \param counter0 pointer to value to update, or 0 to not read
   *  \param counter1 pointer to value to update, or 0 to not read
   *
   *  \returns 0 on success
   *  */ 
  long readCounters(double * counter0, double * counter1 = 0 ); 

  long resetCounters(bool reset0 = true, bool reset1 = false); 

private:
  char* GetName();
  HANDLE _handle; //!
  u3CalibrationInfo _calibration;
  long _error;

  ClassDef(LabJackU3,0)
};


#endif
