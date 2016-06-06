#include "MaxCamPulse.hh"
#include <strings.h>
#include <iostream>
using std::cout;
using std::endl;

ClassImp(MaxCamPulse)
MaxCamPulse::~MaxCamPulse(){}

MaxCamPulse::MaxCamPulse(Int_t nbin)
{
  _nbin=nbin;
  init = true; 
}

void MaxCamPulse::print(ostream & out)
{
  out << "getPulseHeight(): " 
      << getPulseHeight()
      << " getPulseHeightTime(): "
      << getPulseHeightTime()
      << " getPulseStartTime(): "
      << getPulseStartTime()
      << " getPulseStartBin(): "
      << getPulseStartBin()
      << " getPulseEndTime(): "
      << getPulseEndTime()
      << " getPulseEndBin(): "
      << getPulseEndBin()
      << " getPulseIntegral(): "
      << getPulseIntegral()
      << endl;
}
