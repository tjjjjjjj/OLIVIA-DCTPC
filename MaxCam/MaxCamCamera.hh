#ifndef MAXCAM_CAMERA_HH
#define MAXCAM_CAMERA_HH

class TString;
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include <assert.h>

class TH1F;
class TH2S;
class TH2F;
class MaxCamConfig;


class MaxCamCamera  {

public:
    
  // Constructors
  MaxCamCamera();

  virtual ~MaxCamCamera() {};
  
  virtual TString GetName() { return TString("MaxCamCamera"); }
  
  virtual int openCamera(int iCam=0)  { return iCam; }
  virtual int closeCamera() { return -1; }
  virtual int clickCamera() { return -1; }
  virtual int grabImage() { return -1; }
  virtual int expose() { return -1; }
  virtual int cleanCCD()  { return -1; }
  virtual bool imageReady() { return false; }
  
  // manipulate image
  virtual TH2S *createHisto(TString histoName) { histoName+="dummy"; return 0; }
  virtual TH2S *createOverscanHisto(TString histoName) { histoName+="dummy"; return 0; }
  TH2S *getHisto() { return img_histo;}

  // configuration
  MaxCamConfig *config; // Camera configuration.
  MaxCamConfig*  getConfiguration() { return config; }
  virtual int  copyConfiguration() { return -1; }
  virtual int  setHBin(long hbin) {assert(!hbin);return -1; }
  virtual int  setVBin(long vbin) {assert(!vbin);return -1; }
  virtual int  setExposureTime(long exptime);
  virtual int    setTemperature(double temp) { assert(temp>0); return -1; }
  virtual double getTemperature() { return -1; }
  virtual double getGoalTemperature() { return -1; }
  virtual int  setDataBits(int bits) { assert(!bits); return -1; }
  virtual int  getDataBits() {return -1; }
  virtual int  setGain(unsigned short gain) { assert(!gain); return -1; }
  virtual int  getGain() { return -1; }
  virtual int  setDarkFrame() { return -1; }
  virtual int  setNormalFrame() { return -1; }
  virtual int cancelExposure() {return -1; }
  virtual int flushCamera(long repeat=1) { assert(!repeat); return -1; }
  virtual void print() {}
  virtual int  reset(){  return -1; }
  virtual int closeShutter() { return -1; }
  virtual int openShutter() { return -1; }
  
  void deleteBiasFrame();
  int makeBiasFrame(int n=100, bool zeroExposure=false);
  int readBiasFrame(TString fname="");
  void setBiasFrame(TH2F* biasFrame){_biasFrame=biasFrame;}
  void setBiasFrameOverscan(TH2F* biasFrame) {_biasFrameOverscan=biasFrame;}

  TH2S *img_histo; // Histogram of ccd image.
  TH2S *img_histo_overscan; // Histogram of overscan region
  TH2S *img_histo_full; // Histogram of full image (exposed region plus overscan region)

  TH2F *biasFrame() { return _biasFrame; }
  TH2F *_biasFrame; // Bias frame...
  TH2F *biasFrameOverscan() { return _biasFrameOverscan; }
  TH2F *_biasFrameOverscan; // Bias frame...

  int _bnum; // image buffer number  
    
  TStopwatch *_sw; // stopwatch

  virtual int setDigitizeOverscan(bool digi);

  void writeCCDConfigToDB(const char *fname, unsigned int * db_handle = NULL);

private:

  
  int _debugFlag;
  
  ClassDef(MaxCamCamera,0)
    
};

#endif

