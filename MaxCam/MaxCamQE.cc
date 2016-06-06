#include "MaxCamQE.hh"
//#include "MaxCam.hh"

#include <iostream>
#include <fstream>
using std::cout;
using std::endl;
using std::cerr;

// default QE (Kodak KAF-0401E)
float qeff[]={
  250, 0.0,
  260, 0.0,
  270, 0.0,
  280, 0.01,
  290, 0.04,
  300, 0.07,
  310, 0.12,
  320, 0.16,
  330, 0.18,
  340, 0.18,
  350, 0.22,
  360, 0.25,
  370, 0.28,
  380, 0.30,
  390, 0.34,
  400, 0.36,
  410, 0.37,
  420, 0.37,
  430, 0.38,
  440, 0.38,
  450, 0.39,
  460, 0.44,
  470, 0.49,
  480, 0.49,
  490, 0.47,
  500, 0.47,
  510, 0.50,
  520, 0.54,
  530, 0.55,
  540, 0.54,
  550, 0.56,
  560, 0.62,
  570, 0.67,
  580, 0.67,
  590, 0.65,
  600, 0.65,
  610, 0.69,
  620, 0.74,
  630, 0.72,
  640, 0.71,
  650, 0.70,
  660, 0.70,
  670, 0.70,
  680, 0.68,
  690, 0.64,
  700, 0.60,
  710, 0.57,
  720, 0.55,
  730, 0.54,
  740, 0.53,
  750, 0.52,
  760, 0.50,
  770, 0.46,
  780, 0.43,
  790, 0.41,
  800, 0.39,
  810, 0.38,
  820, 0.37,
  830, 0.36,
  840, 0.35,
  850, 0.34,
  860, 0.34,
  870, 0.32,
  880, 0.30,
  890, 0.28,
  900, 0.26,
  910, 0.24,
  920, 0.22,
  930, 0.19,
  940, 0.17,
  950, 0.15,
  960, 0.14,
  970, 0.12,
  980, 0.10,
  990, 0.08,
  1000, 0.07,
  1010, 0.06,
  1020, 0.05,
  1030, 0.04,
  1040, 0.03,
  1050, 0.02,
  1060, 0.02,
  1070, 0.01,
  1080, 0.01,
  1090, 0.01,
  1100, 0.0
};

ClassImp(MaxCamQE)

//____________________
// Class contains quantum efficiency for a CCD camera. Input parameters 
// are a file name that contains columns with wavelengths (nm) and 
// corresponding efficiencies, and a pointer to random number generator.
// If file name is not given, a default QE is used (Kodak KAF-0401E).
// The 2nd argument is needed only in event generators; if not given,
// a local generator will be used (always the same seed).
//
MaxCamQE::MaxCamQE(const char *fname, TRandom *rand) : 
  TGraph(fname) {
  // Constructor uses name of file with QE and a pointer to 
  // the random number generator.
  
  // if random generator not given, use a default
  if (!rand) rnd= new TRandom;
  else rnd = rand;
  
  // if file with QE not given, use a default
  if (!GetN()) {
    int n= sizeof(qeff)/sizeof(float)/2;
    Set( n );
    for (int i=0; i<n; i++) {
      SetPoint(i, qeff[i*2], qeff[i*2+1]);
    }
  }
}


MaxCamQE::MaxCamQE(const MaxCamQE &other) :
  TGraph(other)
{}

MaxCamQE::~MaxCamQE() {
  delete rnd;
}

bool
MaxCamQE::isDetected(float lambda) {

  if (Eval(lambda)>rnd->Rndm()) return true;

  return false;
}


