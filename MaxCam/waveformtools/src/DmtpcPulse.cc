#include "DmtpcPulse.hh"
#include <strings.h>
#include <iostream>
using std::cout;
using std::endl;

ClassImp(DmtpcPulse)

void DmtpcPulse::print(ostream & out)
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
