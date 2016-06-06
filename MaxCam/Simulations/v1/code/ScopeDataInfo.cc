// MaxCam includes
#include "ScopeDataInfo.hh"

// ROOT includes
#include "TObject.h"

// C++ includes
#include <iostream>
using std::cout;
using std::endl;
//________________________________
/*Begin_Html
<center><h2>The SimScopeChannel Class </h2></center>
The SimScopeChannel class is designed to hold the various parameters of the 
Scopes used in the detector.  A separate class is used so that in
future versions of the simulation, several instances of the class may
be added to a TObjectArray or TClonesArray to allow for generation of
multiple Scope signal, each having its own properties.
End_Html*/
//________________________________

ClassImp(ScopeDataInfo)


  ScopeDataInfo::ScopeDataInfo()
{
  //<<<<<Default Constructor>>>>>
  fChanName   = "";
  fChanNumber = 0;
  fConnectedDevice = "";
}
ScopeDataInfo::ScopeDataInfo(Int_t num, TString name, TString cd) 
{
  //<<<<<Alternate Constructor>>>>
  fChanNumber = num;
  fChanName   = name;
  fConnectedDevice = cd;
}

void ScopeDataInfo::setNTriggers(int nt) { _nTriggers = nt; }
int ScopeDataInfo::getNTriggers() const{ return _nTriggers; }

ScopeDataInfo::~ScopeDataInfo()
{
  //<<<<<Destructor>>>>>
}

