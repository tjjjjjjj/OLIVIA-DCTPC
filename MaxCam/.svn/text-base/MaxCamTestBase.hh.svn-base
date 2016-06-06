#ifndef MAXCAM_TEST_BASE_HH
#define MAXCAM_TEST_BASE_HH
#include "TROOT.h"


class MaxCamCamera;
class MaxCamConfig;
class MaxCamChannel;
class TH2F;
class TH1F;
class TDatime;
class TTree;
class TString;
class TFile;
class TCanvas;
class TPad;
class MaxCamTrack;
#include <vector>
using std::vector;
#include <map>
using std::map;
class ScopeHandler;

class MaxCamTestBase  {

public:

  // Constructors
  MaxCamTestBase() {}
  MaxCamTestBase(int debugLevel, TString filePath, TString cameraType);

  virtual ~MaxCamTestBase() {};

  virtual int begin();
  virtual int event();
  virtual int end();

  virtual const char* GetName() { return "MaxCamBase"; }

  void saveEvent();

  MaxCamCamera* ccd(int i=0) { return _camList[i]; }

  virtual bool isTriggered(int icam);
  void setTriggerTrials(int nt) { _triggerTrials=nt; }
  int  getTriggerTrials() { return _triggerTrials; }


  int   makeBiasFrame(int nimages=100);
  int   readBiasFrame(TString fname);
  int   makeBiasFrame(TString fname) { return readBiasFrame(fname); }// historic 
  TH2F* biasFrame(unsigned int icam) { if (icam<0||icam>=_biasFrameList.size()) return 0;   return _biasFrameList[icam]; }
  void deleteBiasFrame();


  TH2F* drawImage(int icam, TString opt="colz", float min=1, float max=0);
  TH1F* drawYields(int icam, TH2F *imhist=0, TString opt="");

  MaxCamChannel* wire() { return _wireHV; }

  MaxCamChannel* mesh() { return _meshHV; }

  MaxCamChannel* pressure() { return _pressure; }

  MaxCamTrack* trfit() { return _trfit; }

  TDatime * dtime() { return _dtime; }

  void saveFile();
  void setSaveFlag(bool doSave=false) { _doSave=doSave; }
  bool getSaveFlag() { return _doSave; }
  const char* getFileName();

  void deleteHotPixels();

  void findHotPixels(int nevents, float th, float rate, long expoTime=200);

  void findHotPixels(TString hotFile);

  void createCanvas();
  TCanvas* getImageCanvas() {return _imageCanvas; }
  TPad* getImagePad() { return _imagePad; }
  TPad* getXProjPad() { return _xprojPad; }
  TPad* getYProjPad() { return _yprojPad; }

  void fillTree(); 

  
  vector<MaxCamCamera*> _camList;
  vector<TH2F*> _biasFrameList;

  ScopeHandler* scope() { return _scopeHandler; }

private:

  ScopeHandler *_scopeHandler;

  TDatime *_dtime; // Time when image was taken.

  TTree *_imageTree; // Tree of events.

  TFile *_file; // Data file

  MaxCamChannel *_wireHV; // Wire HV

  MaxCamChannel *_meshHV; // Mesh HV

  MaxCamChannel *_pressure; // Gas pressure

  MaxCamTrack *_trfit;

  vector<int> _hotPixels;

  double _deltaTemp; // deviation from desired CCD temperature

  int _triggerTrials; // Maximum mumber of images to be tested for trigger.

  bool _doSave; // Automatically save event to file (otherwise use saveEvent method)

  TCanvas *_imageCanvas;
  TPad *_imagePad, *_xprojPad, *_yprojPad;

  ClassDef(MaxCamTestBase,0)
};

#endif

