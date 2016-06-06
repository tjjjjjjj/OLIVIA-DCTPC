#ifndef DMTPC_XCOM_HH
#define DMTPC_XCOM_HH

#include "TGraph.h"

class DmtpcXCOM {

public:

  // Ctors
  //DmtpcXCOM() {};

  DmtpcXCOM(TString filename="");

  //DmtpcXCOM(const DmtpcXCOM &other);
  //virtual ~DmtpcXCOM();

  // x-axis is photon energy in MeV
  // y-axis is attenuation coefficient (mu/rho) in units of cm^2/g
  TGraph* getTotalAttenVsE(bool coherent=true) { 
    return coherent ? _attenVsE_totalWithCoherent : _attenVsE_totalWithoutCoherent;  }

private:

  TGraph* _attenVsE_totalWithCoherent;
  TGraph* _attenVsE_totalWithoutCoherent;
  TGraph* _attenVsE_photoElectric;

  ClassDef(DmtpcXCOM,1)
};

#endif

