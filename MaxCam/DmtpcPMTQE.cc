#include "DmtpcPMTQE.hh"

#include <assert.h>
#include <iostream>
#include <fstream>
using std::cout;
using std::endl;
using std::cerr;

// default QE (just made one up...)
float qe[]={
  200, 0.0,
  250, 0.1,
  260, 0.3,
  280, 0.35,
  300, 0.35,
  500, 0.35,
  700, 0.3,
  800, 0.2,
  900, 0.1,
  1000,0.0
};

ClassImp(DmtpcPMTQE)

//____________________
// Class contains quantum efficiency for a PMT. Input parameters 
// are a file name that contains columns with wavelengths (nm) and 
// corresponding efficiencies, and a pointer to random number generator.
// If file name is not given, a default QE is used.
// If filename = "flat", then a uniform QE will be used, set by val.
// The 3rd argument is needed only in event generators; if not given,
// a local generator will be used (always the same seed).
//
DmtpcPMTQE::DmtpcPMTQE(TString filename, double val, TRandom3 *rand) : 
  TGraph(filename) {
  /* FIXME  if filename is not found then you get a Zombie process */

  // Constructor uses name of file with QE and a pointer to 
  // the random number generator.

  if (filename == "flat") {
    if ((val == -1.0) || (val > 1.0)) {
	assert(!"must provide a value that is <=1.0 to the ctor");
    }
    SetPoint(0, 200, val);   // 200 nm
    SetPoint(1, 1000, val);  // 1000 nm
  }
  
  // if random generator not given, use a default
  if (!rand) rnd= new TRandom3;
  else rnd = rand;
  
  // if file with QE not given, use a default
  if (!GetN()) {
    int n= sizeof(qe)/sizeof(float)/2;
    Set( n );
    for (int i=0; i<n; i++) {
      SetPoint(i, qe[i*2], qe[i*2+1]);
    }
  }
}


DmtpcPMTQE::DmtpcPMTQE(const DmtpcPMTQE &other) :
  TGraph(other)
{}

DmtpcPMTQE::~DmtpcPMTQE() {
  delete rnd;
}

bool
DmtpcPMTQE::isDetected(float lambda) {

  if (Eval(lambda)>rnd->Rndm()) return true;

  return false;
}


