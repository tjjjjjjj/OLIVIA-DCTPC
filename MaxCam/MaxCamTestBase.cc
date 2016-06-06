#include "MaxCamTestBase.hh"
#include "MaxCamConfig.hh"

#ifdef CAM_FLI
#include "MaxCam.hh"
#endif

#ifdef CAM_ALTA
#include "MaxCamAlta.hh"
#endif

#include "MaxCamCamera.hh"
#include "MaxCamChannel.hh"
#include "MaxCamTrack.hh"
#include "MaxCamImageTools.hh"
#include "DmtpcDataConverter.hh"

using namespace DmtpcDataConverter;

//#include "Scope.hh"
//#include "ScopeHandler.hh"

#include <iostream>
#include <fstream>
#include <vector>
using std::cout;
using std::endl;
using std::cerr;
using std::vector;

#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TDatime.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPad.h"

#include<map>
using std::map;
using std::pair;

#include "usb.h"
#define USB_ALTA_VENDOR_ID  0x125c
#define USB_ALTA_PRODUCT_ID  0x0010

ClassImp(MaxCamTestBase)

//____________________
// Class to run MaxCam camera. 
//
//
  MaxCamTestBase::MaxCamTestBase(int debugLevel, TString outputFile, TString cameraType)  {
  // Constructor initializes a camera (one for now) and
  // creates a tree for data-taking.


  cameraType.ToLower();
  bool haveCamera = false;

#ifdef CAM_ALTA
  if (cameraType=="apogee") {
    _camList.push_back( new MaxCamAlta() );
    haveCamera = true;
  }
#endif

#ifdef CAM_FLI
  if (cameraType=="fli") {
    _camList.push_back( new MaxCam(debugLevel) );
    haveCamera = true;
  }
#endif  
  
  assert(haveCamera);
  
  for (unsigned int i=0; i<_camList.size(); i++) {
    assert (!_camList[i]->openCamera(i+1)); 
    // some default setup
    ccd(i)->setExposureTime(0); // exposure
    ccd(i)->setTemperature(-20);   // CCD temperature
  }
  //_cam->setNumberOfFlushes(16); // cleaning of CCD before exposure
  //_cam->startBkgFlush(); // start flushing

  /*
  // initialize scope
  _scopeHandler = new ScopeHandler(Scope::ALAZAR_ATS860);
  scope()->scope()->openScope();
  scope()->scope()->makeScopeConfig(Scope::SCOPE_SETUP_TEST);
  scope()->scope()->initScope();
  */

  _dtime = new TDatime;

  _meshHV   = new MaxCamChannel("_meshHV",    "Mesh Voltage", 0,  0);
  _wireHV   = new MaxCamChannel("_wireHV",    "Wire Voltage", 1,  1);
  _pressure = new MaxCamChannel("_pressure",  "Gas pressure", 2, -1);

  _file = TFile::Open(outputFile,"recreate");
  _imageTree = new TTree("data",  "CCD camera images");
  _imageTree->Branch("timeStamp", "TDatime",       &_dtime,          32000, 0);
  for (unsigned int i=0; i<_camList.size(); i++) {
    TString configName="ccdConfig"; if (i) configName+=i;
    TString imageName="ccdImage"; if (i) imageName+=i;
    _imageTree->Branch( configName, "MaxCamConfig",  &_camList[i]->config,    32000, 0);
    _imageTree->Branch( imageName,  "TH2F",          &_camList[i]->img_histo, 32000, 0);
  }
  _imageTree->Branch("wireHV",    "MaxCamChannel", &_wireHV,   32000, 0);
  _imageTree->Branch("meshHV",    "MaxCamChannel", &_meshHV,   32000, 0);
  _imageTree->Branch("pressure",  "MaxCamChannel", &_pressure, 32000, 0);

  _deltaTemp=3;
  _triggerTrials=10;
  _doSave=false;
}


int
MaxCamTestBase::begin() {

  for (unsigned int i=0; i<_camList.size(); i++) {
    // wait for ccd to cool
    while (ccd(i)->getTemperature()-ccd(i)->getGoalTemperature()>_deltaTemp) {
      cout << GetName() << ": Waiting for CCD " << i << " to cool down, T="<<ccd(i)->getTemperature()
	   <<"  ->  T(goal)=" <<  ccd(i)->getGoalTemperature() << endl;
      gSystem->Sleep(5000);
    }
    cout << ccd(i)->getTemperature() << endl;
    // make one fake readout
    ccd(i)->cleanCCD();
  }

  return 0;
}

int
MaxCamTestBase::readBiasFrame(TString fname) {
  // Read bias frame from a file

  deleteBiasFrame();
  for (unsigned int i=0; i<_camList.size(); i++ ) {
    ccd(i)->clickCamera();
    _biasFrameList.push_back( ccdExpand(ccd(i)->createHisto("bias")) );
    _biasFrameList[i]->Reset();
  }

  TFile biasfile(fname);
  for (unsigned int i=0; i<_camList.size(); i++ ) {
    TString biasFrameName="bias";
    if (i) biasFrameName+=i;
    TH2F* tmpframe=(TH2F*)biasfile.Get(biasFrameName);
    assert (tmpframe);
    _biasFrameList[i]->Add(tmpframe);
  }
  biasfile.Close();

  _file->cd();
  return 0;
}


int
MaxCamTestBase::makeBiasFrame(int nimages) {
  // Take 100 bias frames (no exposure, no shutter) and make one
  // average bias frame. As a result, all subsequent images will have 
  // the bias frame subtracted. In order to stop the bias subtraction
  // call deleteBiasFrame().
  
  cout << endl << GetName() << ": Make bias frame: ";

  deleteBiasFrame();

  for (unsigned int icam=0; icam<_camList.size(); icam++) {
    ccd(icam)->setDarkFrame(); 
    ccd(icam)->setExposureTime(0); 

    for (int i=0; i<nimages; i++) {
      ccd(icam)->clickCamera();
      if (!i) {
	TString biasFrameName="bias";
	if (icam) biasFrameName+=icam;
	_biasFrameList.push_back( ccdExpand(ccd(icam)->createHisto( biasFrameName )) );
      } else { 
	 TH2F *tmpimage = ccdExpand(ccd(icam)->createHisto("addbias"));
	_biasFrameList[icam]->Add( tmpimage );
	delete tmpimage;
      }
    }
    _biasFrameList[icam]->Scale(1./nimages);
  }

  TFile fbias("biasframe.root","RECREATE");
  for (unsigned int icam=0; icam<_camList.size(); icam++) {
    _biasFrameList[icam]->Write();
  }
  fbias.Write();
  fbias.Close();
  
  _file->cd();
  return 0;
}


void
MaxCamTestBase::deleteBiasFrame() {
  // Deletes the bias image so subsequent images are not
  // corrected for the bias.
  
  for (unsigned int icam=0; icam<_biasFrameList.size(); icam++) {
    delete _biasFrameList[icam];
  }
  _biasFrameList.clear();    
}


void
MaxCamTestBase::findHotPixels(TString hotFile) {
  cout << endl << "Find hot pixels: " << hotFile << endl;;
  deleteHotPixels();

  ifstream hf((const char*)hotFile);
  int ipix;
  while (!hf.eof()) {
    hf >> ipix;
    _hotPixels.push_back( ipix );
    cout << ipix << endl;
  }
}


void
MaxCamTestBase::findHotPixels(int nevents, float th, float rate, long expoTime) {


  cout << endl << GetName() << ": Search for hot pixels: ";
  deleteHotPixels();

  ccd()->setNormalFrame();
  ccd()->setExposureTime(expoTime);

  // find pixels over 'th' x RMS in 'nevents'
  map<int,int> hotMap;
  for (int i=0; i<nevents; i++) {

    ccd()->clickCamera();

    TH2F *image  = ccdExpand(ccd()->createHisto("image"));
    TH1F *yields = MaxCamImageTools::createYieldHisto(image);

    float mean = yields->GetMean();
    float rms  = yields->GetRMS();
    int nx = image->GetNbinsX();
    int ny = image->GetNbinsY();
    for (int i=1; i<=nx; i++) {
      for (int j=1; j<=ny; j++) {
	if ( image->GetBinContent(i,j)-mean<rms*th ) continue;
	map<int,int>::iterator ih = hotMap.find(  image->GetBin(i,j) );
	if ( ih==hotMap.end()) hotMap.insert( pair<int,int>( image->GetBin(i,j),1) );
	else ih->second++;
      }
    }
    delete image;
    delete yields;
  }

  // find those repeating and save into ascii file
  ofstream of("hotpixels.dat");
  for (map<int,int>::iterator ih = hotMap.begin(); ih!=hotMap.end(); ih++) {
    if (ih->second<nevents*rate)  continue;
    _hotPixels.push_back( ih->first );
    cout << "hot pixel: " << ih->first << "  rate=" << float(ih->second)/nevents<<endl;
    of << ih->first << endl;
  }
 

}

void
MaxCamTestBase::deleteHotPixels() {
  if (_hotPixels.size()) _hotPixels.clear();
}


TH2F*
MaxCamTestBase::drawImage(int icam, TString opt, float min, float max) {
  // Draw the current image.

  TH2F *tmphist = ccdExpand(ccd(icam)->createHisto("tmphist"));
  if (_biasFrameList[icam]) tmphist->Add( _biasFrameList[icam], -1);
  for (vector<int>::iterator ih=_hotPixels.begin(); ih!=_hotPixels.end(); ih++) {
    tmphist->SetBinContent(*ih,0);
  }
  if (max>min) { tmphist->SetMinimum(min); tmphist->SetMaximum(max); }
  tmphist->DrawCopy(opt);
  return tmphist;
}

TH1F*
MaxCamTestBase::drawYields(int icam, TH2F *imhist, TString opt) {
  // Draw the pixel-yield histogram for the current image.
 
  TH2F* tmphist = imhist ? imhist : drawImage(icam);

  TH1F* hY=MaxCamImageTools::createYieldHisto(tmphist);

  if (opt.Contains("gaus")) hY->Fit("gaus","q");
  return hY;
}


bool
MaxCamTestBase::isTriggered(int icam) {
  // check e.g. that ...

  if (_trfit) delete _trfit;

  TH2F *tmph = ccdExpand(ccd(icam)->createHisto("tmph_tracking"));
  if (biasFrame(icam)) tmph->Add( biasFrame(icam), -1);
  for (vector<int>::iterator ih=_hotPixels.begin(); ih!=_hotPixels.end(); ih++) {
    tmph->SetBinContent(*ih,0);
  }
  TH1F *hY=MaxCamImageTools::createYieldHisto(tmph);
  
  _trfit = new MaxCamTrack(tmph);
  float threshold = hY->GetMean() + hY->GetRMS()*3;
  //cout << "THRRESHOLD="<<threshold<<endl;
  _trfit->setThreshold( threshold);
  _trfit->makeTracks(); 
  cout << GetName() << ": Found # of tracks="<< _trfit->nTracks() << endl;

  //for (int i=0; i<_trfit->nTracks(); i++) {
  //  TH1F *res=_trfit->makeResiduals(i,"y");
  //  cout << "track = " << i << "   residual=" << res->Integral() << endl;
  //  if (res) delete res;
  //  if (_trfit->getTrack(i)) _trfit->getTrack(i)->Draw("same");
  //}
  
  delete hY;

  if ( _trfit->nTracks()<1 || _trfit->nTracks()>3) return false;

  return true;
}



int
MaxCamTestBase::event() {

  for (unsigned int i=0; i<_camList.size(); i++) ccd(i)->getTemperature();
  dtime()->Set();

  int nt=0;
  while(1) {
    
    // this really works up to 2 cameras
    // otherwise, make separate processes
    for (unsigned int i=0; i<_camList.size(); i++) {
      ccd(i)->expose();
    }
    
    /*
    scope()->scope()->acquireTriggers( int(ccd(0)->getConfiguration()->exposureTime*1e-3) );
    if ( !(scope()->scope()->getNValidTriggers() > 0) ) break;
    */

    for (unsigned int i=0; i<_camList.size(); i++) {
      ccd(i)->grabImage ();
    }
    
    if ( ++nt>getTriggerTrials() ) break;

    bool isOk=false;
    for (unsigned int i=0; i<_camList.size(); i++) {
      if (isTriggered(i)) { isOk=true; break; }
    }
    if (isOk) break;
  };

  if (getSaveFlag()) {
    saveEvent();
  }

  return 0;
}

void 
MaxCamTestBase::fillTree() { _imageTree->Fill(); }


void
MaxCamTestBase::saveEvent() {
  // Force saving current event into ROOT file.

  vector<TH2S*> tmplist;
  for (unsigned int i=0; i<_camList.size(); i++) {
    tmplist.push_back( ccd(i)->createHisto("hCCD") );
  }
  fillTree();
  for (unsigned int i=0; i<_camList.size(); i++) {
    delete tmplist[i];
  }
}


int
MaxCamTestBase::end() {
  // Ends datataking and closes root file.

  for (unsigned int i=0; i<_biasFrameList.size(); i++) {
    _biasFrameList[i]->Write();
  }
  saveFile();
  for (unsigned int i=0; i<_camList.size(); i++) ccd(i)->closeCamera();
  cout << GetName() << ": this run has finished and events are saved into " << getFileName() << endl;
  return 0;
}


void
MaxCamTestBase::saveFile() {
  _file->Write(); 
}


const char* 
MaxCamTestBase::getFileName() { return _file->GetName(); }


void
MaxCamTestBase::createCanvas() {
  // Create canvas for image and projections

  _imageCanvas = new TCanvas("imageCanvas","",900,675);
  
  _imagePad = new TPad("imagePad", "", 0.0, 0.2, 0.8, 1.0); _imagePad->SetNumber(0); _imagePad->SetBorderMode(0); _imagePad->Draw();
  _xprojPad = new TPad("xprojPad", "", 0.0, 0.0, 0.8, 0.2); _xprojPad->SetNumber(1); _xprojPad->SetBorderMode(0); _xprojPad->Draw();
  _yprojPad = new TPad("yprojPad", "", 0.8, 0.2, 1.0, 1.0); _yprojPad->SetNumber(2); _yprojPad->SetBorderMode(0); _yprojPad->Draw();
  
  


}
