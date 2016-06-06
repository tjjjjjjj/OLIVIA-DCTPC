//
// LabJackU3.cc
//
#include <iostream>

#include "LabJackU3.hh"
#include "labjackusb.h"

using std::cout;
using std::endl;

ClassImp(LabJackU3)

//LabJackU3::~LabJackU3() {}
LabJackU3::LabJackU3() {
  cout << GetName() << ": constructor" << endl;
  _handle = openUSBConnection(-1);
  //cout << "handle = " << _handle << endl;
  
  //Get calibration information from UE9
  if(getCalibrationInfo(_handle, &_calibration) < 0) {
    if(_error > 0)
      printf("Received an error code of %ld\n", _error);
    closeUSBConnection(_handle);
  }

}

float
LabJackU3::setDAC(Float_t volts, Int_t channel) {
  eDAC(_handle, &_calibration, 0, channel, volts, 0, 0, 0);
  // query and return the actual voltage that was set
  return 0.0;
}

long 
LabJackU3::getUpTime() {
  return getTickCount();
}


long 
LabJackU3::setUpCounters(long pinOffset)
{
  long enableTimers[2], enableCounters[2], timerModes[2]; 
  double timerValues[2];  
  memset(enableTimers,0,sizeof(enableTimers));  
  memset(enableCounters,1,sizeof(enableCounters));  
  memset(timerModes,0,sizeof(timerModes));  
  memset(timerValues,0,sizeof(timerValues));  

  return eTCConfig(_handle, enableTimers, enableCounters, pinOffset, 
                   0,0,timerModes, timerValues, 0, 0); 
}

long 
LabJackU3::readCounters(double * counter0, double * counter1)
{
  long readTimers[2], updateResetTimers[2], readCounters[2], resetCounters[2], ret; 
  double timerValues[2],counterValues[2]; 

  memset(readTimers,0,sizeof(readTimers)); 
  memset(updateResetTimers, 0, sizeof(updateResetTimers)); 
  memset(resetCounters,0,sizeof(resetCounters)); 
        
  readCounters[0] = counter0  ? 1 : 0 ; 
  readCounters[1] = counter1  ? 1 : 0 ; 

  ret = eTCValues (_handle, readTimers, updateResetTimers, readCounters, 
                     resetCounters, timerValues, counterValues, 0, 0) ;  

  if (counter0) *counter0 = counterValues[0]; 
  if (counter1) *counter1 = counterValues[1]; 

  return ret; 
}


long 
LabJackU3::resetCounters(bool reset0, bool reset1)
{
  long readTimers[2], updateResetTimers[2], readCounters[2], resetCounters[2] ; 
  double timerValues[2],counterValues[2]; 

  memset(readTimers,0,sizeof(readTimers)); 
  memset(updateResetTimers, 0, sizeof(updateResetTimers)); 
  memset(readCounters,0,sizeof(readCounters)); 

  resetCounters[0] = reset0; 
  resetCounters[1] = reset1; 
        
  return eTCValues (_handle, readTimers, updateResetTimers, readCounters, 
                     resetCounters, timerValues, counterValues, 0, 0) ;  
}


char* 
LabJackU3::GetName() { return "LabJackU3"; }


