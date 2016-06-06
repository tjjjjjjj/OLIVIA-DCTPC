#ifndef MAXCAM_DATASET_HH
#define MAXCAM_DATASET_HH

#include <assert.h>

#include "TROOT.h"
#include "TFile.h"

class MaxCamConfig;
class MaxCamChannel;
class TChain;
class TTree;
class TDatime;
class TH2F;

class MaxCamDataset  {

public:

  // Ctors
    MaxCamDataset( const char* fileName="dataset.root", TString foption="");

    MaxCamDataset(const MaxCamDataset &other);

    virtual ~MaxCamDataset();

    TString GetName() { return "MaxCamDataset"; }
    
    void Write() { _file->Write(); }

    MaxCamChannel *wire() { return _wireHV; }
    
    MaxCamChannel *mesh() { return _meshHV; }

    MaxCamChannel *pressure() { return _pressure; }

    TTree* tree() { return _imageTree; }

    TChain* chain() { return (TChain*)_imageTree; }

    MaxCamConfig* ccdConfig() { return _config; }

    TDatime* timeStamp() { return _dtime; }

    TH2F* ccdImage() { return _img_histo; }

    void fill(); 

    
//private:
    
    TTree* _imageTree; // tree with events
    TFile* _file; // event file
    TH2F*  _img_histo; // image histogram
    MaxCamChannel* _meshHV; // dift voltage
    MaxCamChannel* _wireHV; // anode voltage
    MaxCamChannel* _pressure; // pressure
    MaxCamConfig * _config; // ccd configuration
    TDatime *_dtime; // time-stamp
    
    ClassDef(MaxCamDataset,0)

        };

#endif

