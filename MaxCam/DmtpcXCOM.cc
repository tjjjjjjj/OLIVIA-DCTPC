#include "DmtpcXCOM.hh"

#include <assert.h>
#include <iostream>
#include <fstream>
using std::cout;
using std::endl;
using std::cerr;

ClassImp(DmtpcXCOM)

//____________________
// Class contains 
//
DmtpcXCOM::DmtpcXCOM(TString filename) {

  TString fmtTotalWithCoherent = "%lg %*lg %*lg %*lg %*lg %*lg %lg %*lg";
  TString fmtTotalWithoutCoherent = "%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg";
  TString fmtPhotoelectric = "%lg %*lg %*lg %lg";

  _attenVsE_totalWithCoherent = new TGraph(filename, fmtTotalWithCoherent);
  _attenVsE_totalWithoutCoherent = new TGraph(filename, fmtTotalWithoutCoherent);
  _attenVsE_photoElectric = new TGraph(filename, fmtPhotoelectric);

}


//DmtpcXCOM::DmtpcXCOM(const DmtpcXCOM &other) :
//  TGraph(other)
//{}

//DmtpcXCOM::~DmtpcXCOM() {
//  delete rnd;
//}

