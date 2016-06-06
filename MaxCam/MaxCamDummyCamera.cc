#include "MaxCamDummyCamera.hh"
#include "MaxCamConfig.hh"
#include <iostream>
#include <fstream>
#include <vector>
using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::flush;

#include "TString.h"
#include "TFile.h"
#include "TH2S.h"
#include "TStopwatch.h"



ClassImp(MaxCamDummyCamera)

//____________________
// Class to run Dummy Camera. 
//
// *** Example:
//
//
MaxCamDummyCamera::MaxCamDummyCamera(int debugLevel)  {
  // Constructor

  config = new MaxCamConfig("ccdconfig","CCD Camera configuration");;
  img_histo = new TH2S;
  
  _sw = new TStopwatch;
}



int 
MaxCamDummyCamera::openCamera(int i) {
  // Prepare camera for taking data: estalish a connection and set the image area;

  cout << GetName() << ": configuration of dummy camera " << i << endl;
  return 0;
}



TH2S*
MaxCamDummyCamera::createHisto(TString histoName) {
  // Create a histrogram from the image array in the memory.

  int nx = 1;
  int ny = 1;

  img_histo = new TH2S(histoName, "", 
		       nx, 0.0, 1024. ,  
		       ny, 0.0, 1024. );

  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      img_histo->SetBinContent(i, j, 1);
    }
  }

  return img_histo;
}
