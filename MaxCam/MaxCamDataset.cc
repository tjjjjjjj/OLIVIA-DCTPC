#include "MaxCamDataset.hh"
#include "MaxCamChannel.hh"
#include "TDatime.h"
#include "MaxCamConfig.hh"
#include "TChain.h"
#include "TMath.h"
#include "TH2F.h"

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;




ClassImp(MaxCamDataset)

//____________________
//
    MaxCamDataset::MaxCamDataset( const char *fname, TString foption)  {
    // 
    
    _dtime = new TDatime;

    _meshHV   = new MaxCamChannel("_meshHV",    "Mesh Voltage", 0, -1);
    _wireHV   = new MaxCamChannel("_wireHV",    "Wire Voltage", 1, -1);
    _pressure = new MaxCamChannel("_pressure",  "Gas pressure", 2, -1);
    
    _config=new MaxCamConfig;
    
    _img_histo=new TH2F("image","",96,0,768, 64,0,512);
    
    foption.ToLower();
    
    if (fname && foption=="recreate") {
        _file = TFile::Open(fname, foption);
        _imageTree = new TTree("data",  "CCD camera images");
        tree()->Branch("timeStamp", "TDatime",       &_dtime,     32000, 0);
        tree()->Branch("ccdConfig", "MaxCamConfig",  &_config,    32000, 0);
        tree()->Branch("ccdImage",  "TH2F",          &_img_histo, 32000, 0);
        tree()->Branch("wireHV",    "MaxCamChannel", &_wireHV,    32000, 0);
        tree()->Branch("meshHV",    "MaxCamChannel", &_meshHV,    32000, 0);
        tree()->Branch("pressure",  "MaxCamChannel", &_pressure,  32000, 0);
    }
    else if (fname && foption=="") {
        _file = TFile::Open(fname, foption);
        _imageTree = new TChain("data",  "CCD camera images");
        chain()->Add(fname);
        chain()->SetBranchAddress("timeStamp", &_dtime);
        chain()->SetBranchAddress("ccdConfig", &_config);
        chain()->SetBranchAddress("ccdImage",  &_img_histo);
        chain()->SetBranchAddress("wireHV", &_wireHV);
        chain()->SetBranchAddress("meshHV", &_meshHV);
        chain()->SetBranchAddress("pressure", &_pressure);            
    }
    else if (fname) assert(0);
    
}


MaxCamDataset::MaxCamDataset(const MaxCamDataset &other) {
    cout << GetName() << "Copy Constructor not done" << endl;
    _file = other._file;
}

MaxCamDataset::~MaxCamDataset() {
    cout << GetName() << "Destructor not done" << endl;
}


void
MaxCamDataset::fill() { _imageTree->Fill(); }
