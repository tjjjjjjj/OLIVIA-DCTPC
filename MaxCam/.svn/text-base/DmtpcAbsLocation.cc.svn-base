#include "DmtpcAbsLocation.hh"
#include "TROOT.h"
#include "TVector3.h"
#include <fstream>
#include <iostream>
using namespace std;

ClassImp(DmtpcAbsLocation)

DmtpcAbsLocation::DmtpcAbsLocation(double lat, double lon, double z, double alpha, double beta, double gamma) {

   _position = new TVector3(lat,lon,z);
   _angles = new TVector3(alpha,beta,gamma);
  
}

DmtpcAbsLocation::DmtpcAbsLocation(TVector3* position, TVector3* angles) {
   
   _position = new TVector3(position->X(),position->Y(),position->Z());
   _angles = new TVector3(angles->X(),angles->Y(),angles->Z());
   
}

DmtpcAbsLocation::~DmtpcAbsLocation() {
  cout << GetName() << "Destructor not done" << endl;
}

