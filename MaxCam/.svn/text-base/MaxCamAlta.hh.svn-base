#ifndef MAXCAMALTA_HH
#define MAXCAMALTA_HH

#include "TROOT.h"
#include "MaxCamCamera.hh"

class TString;
class TH1F;
class TH2S;
class CApnCamera;


class MaxCamAlta : public MaxCamCamera {

#define N_OVERSCAN_COLS 8
  //#define N_OVERSCAN_ROWS 0  // have not implemented this yet

public:

  // Constructors
  MaxCamAlta(int debug=0);

  virtual ~MaxCamAlta() {};

  virtual TString GetName() { return TString("MaxCamAlta"); }

  // Camera action

  virtual int openCamera(int iCam=0);
  virtual int closeCamera();
  virtual int clickCamera();
  virtual int grabImage();
  virtual int expose();
  virtual int cleanCCD() { return 0; }
  virtual int flushCamera(long repeat=1);
  virtual bool imageReady();
  virtual int closeShutter();
  virtual int openShutter();

  virtual void print();
  virtual int  reset();
  int resetWithFlush();

  // manipulate image
  virtual TH2S *createHisto(TString histoName);
  virtual TH2S *createOverscanHisto(TString histoName);
  virtual TH2S *createFullHisto(TString histoName);

  // configuration
  virtual int copyConfiguration();
  virtual int  setHBin(long hbin);
  virtual int  setVBin(long vbin);
  virtual int  setExposureTime(long exptime);
  virtual int  setTemperature(double temp);
  virtual double getTemperature();
  virtual double getGoalTemperature();
  virtual int  setDataBits(int bits);
  virtual int  getDataBits();
  virtual int  setGain(unsigned short gain);
  virtual int  getGain();
  virtual int  setDarkFrame();
  virtual int  setNormalFrame();
  int setTestFrame();
  virtual int  cancelExposure();
  int setFanSpeed(int speed);

  static void discoverCameras();
  static unsigned short  nCamera;
  
  virtual int setDigitizeOverscan(bool digi);

  CApnCamera* getAlta() { return _alta; }

private:

  CApnCamera *_alta; // Pointer to camera driver


  ClassDef(MaxCamAlta,0)
};

#endif
