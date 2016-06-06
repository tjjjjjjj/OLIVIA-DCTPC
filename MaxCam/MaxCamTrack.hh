#ifndef MAXCAM_TRACK_HH
#define MAXCAM_TRACK_HH

#include "TObject.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TVector2.h"
#include <vector>
using std::vector;

class TCanvas;

class MaxCamTrack : public TObject {

public:

  // Ctors

    MaxCamTrack();

    MaxCamTrack(TH2F *image, bool rotate90deg=false);

    MaxCamTrack(const MaxCamTrack &other);

  virtual ~MaxCamTrack();

  void setImage( TH2F* image) { _image = image; }

  void setThreshold(float th) { _th = th; }

  void checkThresholdAndNeighbors();
  void findMaxBin();
  void findSlope(int i=-1, int j=-1);
  void fitTrack();
  void cleanUpPixels();

  TF1 *getTrack() { return _ftrack; }

  TH2F *getImage() { return _image; }
  TH2F *getFitImage() { return _imageth; }

  vector<TF1*> trackList;
  int nTracks() { return trackList.size(); }
  TF1* getTrack(int i) { if (i<0 || i>=nTracks()) return 0;  return trackList[i]; }

  void makeTracks(TCanvas *dbgCanvas=0);
  TH1F* makeResiduals(int itr, TString opt, float width=-1);

  TH1F* getAngleHisto() { return _hphi; }

  void setWireBinList( vector<int> wbl, TString wc) { _wireBinList=wbl; setWireCoordinate(wc);  }
  vector<int> getWireBinList() { return _wireBinList; }
  
  void setWireCoordinate(TString wc) { _wireCoord=wc; }
  TString getWireCoordinate() { return _wireCoord; }
  
  void setHoughHistogram(int nb, float min, float max) {
    _hphibins=nb;
    _hphimin=min;
    _hphimax=max;
  }

  void plotDebug(TCanvas *dbgCanvas);

  void clearMemory();

  static TVector2 distanceFromPoint(TF1 *t, double x, double y);

  void setCleaningSpan(float dist) { _yroi=dist; }

private:

    TH1F *_hphi; // phi histogram
    
    TH2F *_image; // Pointer to original (calibrated) image
    
    TH2F *_imageth; // Image after applied threshold.
    
    TH2F *_imagephi; // Array of phi values for each pixel
    
    float _th; // Threshold for accepting pixel in track finding
    
    int   _npixel, _nx, _ny; 
    float _xratio; 
    float _yratio;
    
    
    TF1 *_ftrack; // Fitted track
    
    int _i0, _j0;
    
    vector<int> _wireBinList;

    vector<int> _peakList;
    bool isNew(int ix, int iy);
    
    TString _wireCoord; // wire coordinate (x or y - assumes wires are parallel either with y or x axis)

    int _hphibins; // number of bins for Hough transform histogram
    float _hphimin, _hphimax; // min, max for Hough transform histogram
    double _pihalf;
    float _yroi;  // Size of y band for selecting pixels.
    int _phiroi; // Width of phi histogram used to compute mean slope at track search stage.
    
    ClassDef(MaxCamTrack,1)
};

#endif

