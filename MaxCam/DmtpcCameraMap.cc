#include "DmtpcCameraMap.hh"
#include "DmtpcStringTools.hh"

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "TMath.h"
using std::cout;
using std::endl;
using std::cerr;
using std::ifstream;
using std::string;
using std::getline;
using std::stringstream;

ClassImp(DmtpcCameraMap)

DmtpcCameraMap::DmtpcCameraMap()  {
  _verbose = 0;
}

DmtpcCameraMap::~DmtpcCameraMap() {
}

void 
DmtpcCameraMap::addCamera(TString serialNumber, Int_t pad, TString rotationType) {
  _serialNumbers.push_back(serialNumber);
  _pads.push_back(pad);
  _rotations.push_back(rotationType);
  return;
}


void 
DmtpcCameraMap::loadMap(TString filename) {

  TString snX, rotX;
  Int_t padX;

  string line;
  ifstream configfile(filename.Data());
  if (configfile.is_open()) {
    while ( !configfile.eof() ) {
      getline(configfile, line);
      //configfile >> snX >> padX >> rotX;

      // trim leading and trailing whitespace
      DmtpcStringTools::trim(line);

      // ignore empty lines
      if (line.size() == 0) continue;

      // ignore comment lines (lines whose first non-whitespace char is #)
      if (line.find("#") == 0) continue;

      if (_verbose > 0) cout << line << endl;
      stringstream ss;
      ss << line;
      ss >> snX >> padX >> rotX;
      if (_verbose > 0) {
	cout << "SN:   " << snX << endl;
	cout << "Pad:  " << padX << endl;
	cout << "Rot:  " << rotX << endl;
	cout << endl;
      }
      addCamera(snX, padX, rotX);
    }
  }
  configfile.close();
  return;
}

Int_t 
DmtpcCameraMap::getPad(TString sn) {
  return _pads.at(getIdOfSN(sn));
}

TString
DmtpcCameraMap::getRot(TString sn) {
  return _rotations.at(getIdOfSN(sn));
}

double DmtpcCameraMap::getRotRadians(TString sn)
{

  if (getRot(sn) == "QUARTER_LEFT") return -TMath::Pi()/2; 
  if (getRot(sn) == "QUARTER_RIGHT") return TMath::Pi()/2; 
  if (getRot(sn) == "HALF_TURN") return TMath::Pi(); 

  return 0; 
}

Int_t
DmtpcCameraMap::getIdOfSN(TString sn) {
  Int_t nelements = _serialNumbers.size();
  
  for (int ii=0; ii<nelements; ii++) {
    if (_serialNumbers.at(ii) == sn) {
      return ii;
    }
  }
  cout << "DmtpcCameraMap.getIdOfSN:  no match found... returning -1" << endl;
  return -1;

}
