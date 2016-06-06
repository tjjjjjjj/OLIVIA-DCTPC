#ifndef MAXCAM_HH
#define MAXCAM_HH

#include "TROOT.h"
#include "MaxCamCamera.hh"

typedef long int fliframe_t;
typedef long int flidev_t;
typedef unsigned short u_int16_t;
class TString;

 

class MaxCam : public MaxCamCamera {

public:

  // Constructors
  MaxCam(int debugLevel=0);

  virtual ~MaxCam() {};

  // Camera action

  virtual int openCamera(int iCam=0);
  virtual int closeCamera();
  virtual int clickCamera();
  virtual int expose();
  virtual int grabImage();

  // Save current image in various formats
  int writePNG(TString fName);
  int writeFITS(TString fName);
  int writeROOT(TString fName);

  // manipulate image
  virtual TH2S *createHisto(TString histoName);
  TH1F *createYieldHisto(TString histoName);


  // Camera info

  void  printCameraInfo(int iCam=0);
  char* getCameraModel();
  long  getHWRevision();
  long  getFWRevision();  

  int numcams; // number of cameras attached to the system

  // Configure camera

  virtual int setExposureTime(long exptime);

  long getArrayArea  (TString what="ULX");
  long getVisibleArea(TString what="ULX");
  int  setImageArea(long ul_x=-1, long ul_y=-1, long lr_x=-1, long lr_y=-1);
  virtual int  setHBin(long hbin);
  virtual int  setVBin(long vbin);

  virtual int    setTemperature(double temp);
  virtual double getTemperature();
  virtual double getGoalTemperature();
  //double readTemperature();

  virtual int cleanCCD();
  int setNumberOfFlushes(long nflush);
  int startBkgFlush();
  int stopBkgFlush();
  virtual int flushCamera(long repeat=1);

  virtual int closeShutter();
  virtual int openShutter();
  int triggerShutter();

  int totalExposures;

  virtual int setDarkFrame();
  virtual int setNormalFrame();

  double getPixelSize(TString what="X");

  virtual int cancelExposure();


private:

  u_int16_t *img; // Pointer to camera image.

  TStopwatch *_sw;

  void findCams(); 

  ClassDef(MaxCam,0)
};

#endif


/*
LIBFLIAPI FLIOpen(flidev_t *dev, char *name, flidomain_t domain);
LIBFLIAPI FLISetDebugLevel(char *host, flidebug_t level);
LIBFLIAPI FLIClose(flidev_t dev);
LIBFLIAPI FLIGetLibVersion(char* ver, size_t len);
LIBFLIAPI FLIGetModel(flidev_t dev, char* model, size_t len);
LIBFLIAPI FLIGetPixelSize(flidev_t dev, double *pixel_x, double *pixel_y);
LIBFLIAPI FLIGetHWRevision(flidev_t dev, long *hwrev);
LIBFLIAPI FLIGetFWRevision(flidev_t dev, long *fwrev);
LIBFLIAPI FLIGetArrayArea(flidev_t dev, long* ul_x, long* ul_y,
			  long* lr_x, long* lr_y);
LIBFLIAPI FLIGetVisibleArea(flidev_t dev, long* ul_x, long* ul_y,
			    long* lr_x, long* lr_y);
LIBFLIAPI FLISetExposureTime(flidev_t dev, long exptime);
LIBFLIAPI FLISetImageArea(flidev_t dev, long ul_x, long ul_y,
			  long lr_x, long lr_y);
LIBFLIAPI FLISetHBin(flidev_t dev, long hbin);
LIBFLIAPI FLISetVBin(flidev_t dev, long vbin);
LIBFLIAPI FLISetFrameType(flidev_t dev, fliframe_t frametype);
LIBFLIAPI FLICancelExposure(flidev_t dev);
LIBFLIAPI FLIGetExposureStatus(flidev_t dev, long *timeleft);
LIBFLIAPI FLISetTemperature(flidev_t dev, double temperature);
LIBFLIAPI FLIGetTemperature(flidev_t dev, double *temperature);
LIBFLIAPI FLIGrabRow(flidev_t dev, void *buff, size_t width);
LIBFLIAPI FLIExposeFrame(flidev_t dev);
LIBFLIAPI FLIFlushRow(flidev_t dev, long rows, long repeat);
LIBFLIAPI FLISetNFlushes(flidev_t dev, long nflushes);
LIBFLIAPI FLISetBitDepth(flidev_t dev, flibitdepth_t bitdepth);
LIBFLIAPI FLIReadIOPort(flidev_t dev, long *ioportset);
LIBFLIAPI FLIWriteIOPort(flidev_t dev, long ioportset);
LIBFLIAPI FLIConfigureIOPort(flidev_t dev, long ioportset);
LIBFLIAPI FLILockDevice(flidev_t dev);
LIBFLIAPI FLIUnlockDevice(flidev_t dev);
LIBFLIAPI FLIControlShutter(flidev_t dev, flishutter_t shutter);
LIBFLIAPI FLIControlBackgroundFlush(flidev_t dev, flibgflush_t bgflush);
LIBFLIAPI FLISetDAC(flidev_t dev, unsigned long dacset);
LIBFLIAPI FLIList(flidomain_t domain, char ***names);
LIBFLIAPI FLIFreeList(char **names);
LIBFLIAPI FLISetFilterPos(flidev_t dev, long filter);
LIBFLIAPI FLIGetFilterPos(flidev_t dev, long *filter);
LIBFLIAPI FLIGetFilterCount(flidev_t dev, long *filter);
LIBFLIAPI FLIStepMotor(flidev_t dev, long steps);
LIBFLIAPI FLIStepMotorAsync(flidev_t dev, long steps);
LIBFLIAPI FLIGetStepperPosition(flidev_t dev, long *position);
LIBFLIAPI FLIGetStepsRemaining(flidev_t dev, long *steps);
LIBFLIAPI FLIHomeFocuser(flidev_t dev);
LIBFLIAPI FLICreateList(flidomain_t domain);
LIBFLIAPI FLIDeleteList(void);
LIBFLIAPI FLIListFirst(flidomain_t *domain, char *filename,
		      size_t fnlen, char *name, size_t namelen);
LIBFLIAPI FLIListNext(flidomain_t *domain, char *filename,
		      size_t fnlen, char *name, size_t namelen);
LIBFLIAPI FLIReadTemperature(flidev_t dev,
					flichannel_t channel, double *temperature);
LIBFLIAPI FLIGetFocuserExtent(flidev_t dev, long *extent);
*/
