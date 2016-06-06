#ifndef MAXCAM_DUMMY_CAMERA_HH
#define MAXCAM_DUMMY_CAMERA_HH

#include "TROOT.h"
#include "MaxCamCamera.hh"

class TString;
class TH1F;
class TH2S;

class MaxCamDummyCamera : public MaxCamCamera {

public:

  // Constructors
  MaxCamDummyCamera(int debugLevel=0);

  virtual ~MaxCamDummyCamera() {};

  virtual int openCamera(int iCam=0);
  virtual TString GetName() { return TString("MaxCamDummyCamera"); }

  // Camera action

  // manipulate image
  virtual TH2S *createHisto(TString histoName);

  virtual bool imageReady() { return true; }

private:

  ClassDef(MaxCamDummyCamera,0)
};

#endif
