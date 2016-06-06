
#ifndef MAXCAMANDOR_HH
#define MAXCAMANDOR_HH

#include "TROOT.h"
#include "MaxCamCamera.hh"


class TString;
class TH1F;
class TH2F;
class CApnCamera;


class MaxCamAndor : public MaxCamCamera {

#define N_OVERSCAN_COLS 8
  //#define N_OVERSCAN_ROWS 0  // have not implemented this yet

public:

  // Constructors
  MaxCamAndor(int debug=0);

  virtual ~MaxCamAndor() {};

  virtual TString GetName() { return TString("MaxCamAndor"); }

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

  // Andor specific
  int setFanSpeed(int speed);
  int coolerOn();
  int coolerOff();

  int baselineClampOn();
  int baselineClampOff();
  
  float getReadoutTime();
  
  void setImage();
  void calcDaqTime();

  int setPreAmpGain(int i);
  
  virtual int setDigitizeOverscan(bool digi);


private:
  long *imageData;

  ClassDef(MaxCamAndor,0)
};

#endif
