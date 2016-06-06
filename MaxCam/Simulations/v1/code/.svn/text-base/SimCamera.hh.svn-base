#ifndef __SIMCAMERA__
#define __SIMCAMERA__

#include "TH2.h"
#include "TObject.h"
#include "TRandom3.h"
#include "TCutG.h"

using namespace::std;

class SimCamera : public TObject {

public: 

  SimCamera();
  SimCamera(double x, double y, int xbin, int ybin,int pixperbin,
	    double xwidth,double ywidth,int camnum); 
  //SimCamera(const SimCamera &other);
  virtual ~SimCamera();

  void setSerialNumber(TString ser){fSerialNumber = ser;}
  void setCameraNumber(int num){fCameraNumber = num;}
  void setPosition(double x, double y);
  void setWidths(double xwidth,double ywidth);
  void setBins(int xbin, int ybin);
  void setPixelPerBin(int pixperbin){fPixelsPerBin = pixperbin;}
  void setEMGain(double emgain){fEmGain = emgain;}
  void setNoiseFactor(double noisefactor){fNoiseFactor = noisefactor;}
  void setReadNoise(double noise){fReadNoise = noise;}
  void setDarkCurrent(double current){fDarkCurrent = current;}
  void setGain(double gain){fGain = gain;}
  void setGain_Landau(double landau_1,double landau_2);
  void setRadialParams(double radial1,double radial2,double radial3);
  void setPhotonsADU(double photonsadu);
  void setBias(double bias){fBias = bias;}
  void setGainMap(TH2F* gainmap) { fGainMap = (TH2F*) gainmap->Clone();} 
  void setActiveRegion(TCutG *active) { fActiveRegion = (TCutG*) active->Clone();} 
  void setActiveRegionUnits(TString actunits) { fActiveRegionUnits=actunits; }
  void setCountBeforeSpacers(double count){fCountBfSpacers = count;}
  void setCountBeforeRadialEffect(double count){fCountBfRadial = count;}
  void setFinalCount(double count){fFinalCount = count;}
  void setImage();
  void resetImage();
  void applyNoise();
  void emptyImage();

  void normalizeGainMap(); 

  TString getSerialNumber(){return fSerialNumber;}
  int getCameraNumber(){return fCameraNumber;}
  double getPositionX(){return fX;}
  double getPositionY(){return fY;}
  double getWidthX(){return fXWidth;}
  double getWidthY(){return fYWidth;}
  int getBinsX(){return fXBins;}
  int getBinsY(){return fYBins;}
  int getPixelsX(){return fXBins*fPixelsPerBin;}
  int getPixelsY(){return fYBins*fPixelsPerBin;}
  int getPixelPerBin(){return fPixelsPerBin;}
  double getEMGain(){return fEmGain;}
  double getNoiseFactor(){return fNoiseFactor;}
  double getReadNoise(){return fReadNoise;}
  double getDarkCurrent(){return fDarkCurrent;}
  double getGain(){return fGain;}
  double getGain(int u, int v);
  double getGainLandau_1(){return fGainLandau_1;}
  double getGainLandau_2(){return fGainLandau_2;}
  double getRadialParams_1(){return fRadialParams_1;}
  double getRadialParams_2(){return fRadialParams_2;}
  double getRadialParams_3(){return fRadialParams_3;}
  double getPhotonsADU(){return fPhotonsADU;}
  double getBias(){return fBias;}
  double getBinPositionX(int i, TString opt = "");
  double getBinPositionY(int i, TString opt = "");
  double getPixelPositionX(int i, TString opt = "");
  double getPixelPositionY(int i, TString opt = "");
  double getPixelX(double x);
  double getPixelY(double y);
  int getBinX(double x){return static_cast<int>(getPixelX(x))/fPixelsPerBin+1;}
  int getBinY(double y){return static_cast<int>(getPixelY(y))/fPixelsPerBin+1;}
  double getCountBeforeSpacers(){return fCountBfSpacers;}
  double getCountBeforeRadialEffect(){return fCountBfRadial;}
  double getFinalCount(){return fFinalCount;}
  bool isInImage(double x,double y);
  TH2F* getCCDImage(){return fccdImage;}
  TH2S* getRawCCDImage();

private:

  TString fSerialNumber; //Camera Serial Number
  int fCameraNumber;     //Camera Number
  double fX;             //x-Position of camera center in mm
  double fY;             //y-Position of camera center in mm
  int fXBins;            //Number of bins on x-axis
  int fYBins;            //Number of bins on y-axis
  int fPixelsPerBin;     //Width of each bin in pixel widths
  double fXWidth;        //x Width of image in mm
  double fYWidth;        //y Width of image in mm
  double fEmGain;        //EM gain of EMCCD camera
  double fNoiseFactor;   //Extra noise factor of EMCCD camera
  double fReadNoise;     //Gaussian read noise in count/bin
  double fDarkCurrent;   //Gaussian dark current noise in count/bin
  double fGain;          //Total gain at center of image in count/keV
  double fGainLandau_1;  //gain smearing via Landau distribution, first parameter
  double fGainLandau_2;  //gain smearing via Landau distribution, second parameter
  double fRadialParams_1; 
  double fRadialParams_2;
  double fRadialParams_3;
  double fPhotonsADU;
  double fBias;          //Bias in counts
  double fCountBfSpacers;//Count before spacer/radial effects
  double fCountBfRadial; //Count after spacers, before radial effect
  double fFinalCount;    //Count after spacers, radial effect
  TH2F *fccdImage;       //Image histogram
  TH2F *fGainMap;        //Gain Map 
  TRandom3 *fRnd;        //Random number generator
  TCutG *fActiveRegion;  //Active region 
  TString fActiveRegionUnits; // units of active region points

  ClassDef(SimCamera,3)
};

#endif
