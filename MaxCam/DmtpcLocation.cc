#include "DmtpcLocation.hh"
#include "TROOT.h"
#include "TVector3.h"
#include <fstream>
#include <iostream>
using namespace std;

ClassImp(DmtpcLocation)

DmtpcLocation::DmtpcLocation(double x, double y, double z, double alpha, double beta, double gamma) {

   _position = new TVector3(x,y,z);
   _angles = new TVector3(alpha,beta,gamma);
  
}

DmtpcLocation::DmtpcLocation(TVector3* position, TVector3* angles) {
   
   _position = new TVector3(position->X(),position->Y(),position->Z());
   _angles = new TVector3(angles->X(),angles->Y(),angles->Z());
   
}

DmtpcLocation::~DmtpcLocation() {
  cout << GetName() << "Destructor not done" << endl;
}

