#ifndef DMTPC_EVENT_VIEWER_HH
#define DMTPC_EVENT_VIEWER_HH

#include "TROOT.h"
#include "TString.h"

class TCanvas;
class TH2F;
class DmtpcDataset;
class DmtpcCameraMap;
class TTree;

class DmtpcEventViewer {
    
public:
    
  DmtpcEventViewer();
  DmtpcEventViewer(TString fname, TString whatToPlot="scope ccd", Int_t maxccd=2);
  virtual ~DmtpcEventViewer();
  
  void next();
  void prev();
  void validateEventNumber();
  void dark();
  void step(int pause_ms=100, int evMin=-1, int evMax=-1);

  TTree *tree();
  
  void show(int ii);
  
  void getEvent(int iev);
  
  DmtpcDataset *data();
  
  void setRange(float min, float max) {
    minRange=min;
    maxRange=max;
  }
  void setSigma(float sig) { _nsigma = sig; }

  void setDarkSubtract(bool ds) { _darkSubtract = ds; }

  Int_t getMaxCCD() { return maxCCD; }
  void setVerbose(Int_t verb) { verbose = verb; }
  Int_t getVerbose() { return verbose; }
  
  void setNScope(int n) { _nscope=n; }
  
  void loadCameraLayoutMap();
  DmtpcCameraMap* cameraMap() { return _cm; }
  Int_t getPadNumber(Int_t iccd);
  
  void drawAnodeBoundary();

  TH2F* img(int iccd) { return image[iccd]; }
  TString sn(int iccd);

protected:
  void display();
  void displayCCD();
  void displayScope();
  


private:
  
  DmtpcDataset *d; // DMTPC dataset
  
  TH2F *image[4], *bias[4];
  
  TString what; // what to display
  
  int iev; // current event number
  
  TCanvas *c_ccd; // camvas for displaying ccd images
  TCanvas *c_scope; // canvas for displaying waveforms
  
  
  int maxWF; // maximum nummber of waveforms to display
  
  int maxCCD; // maximum number of ccd images to display
  
  float minRange; // minimum yield to display
  float maxRange; // maximum yield to display
  
  Int_t verbose;
  
  int _nscope; // number of scope boards
  
  float _nsigma;
  bool _darkSubtract;
  TString _cameraMapName;
  bool _cameraMapExists;
  DmtpcCameraMap *_cm;

  ClassDef(DmtpcEventViewer,0)
    
};


#endif
