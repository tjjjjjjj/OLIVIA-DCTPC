#include "DmtpcEventViewer.hh"
#include "DmtpcCameraMap.hh"
#include "DmtpcDataset.hh"
#include "MaxCamConfig.hh"
#include "MaxCamWaveformTools.hh"
#include "MaxCamImageTools.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "ScopeDataInfo.hh"
#include "TCanvas.h"
#include "McDarkTrack.hh"
#include "TStyle.h"
#include "TSystem.h"
#include "TArc.h"

#include <iostream>
using std::cout;
using std::endl;

ClassImp(DmtpcEventViewer)


//____________________
// A simple event viewer for the DMTPC data. 
//
// DmtpcEventViewer ev("file_name.root");
// ev.next();
// ev.next();
// ...
//
//
//    
  DmtpcEventViewer::DmtpcEventViewer(TString fname, TString whatToPlot, Int_t maxccd) :
    what(whatToPlot),
    iev(-1),
    maxWF(8),
    maxCCD(maxccd),
    minRange(-50),
    maxRange(250),
    _nscope(2),
    _nsigma(3.0),
    _darkSubtract(true),
    _cameraMapName("camera.map"),
    _cameraMapExists(false) {

    gStyle->SetOptTitle(1);
    
    if (what.Contains("ccd")) { 
      c_ccd= new TCanvas("c_ccd","c_ccd", 600, 400); 
      if (maxCCD != 1) c_ccd->Divide(2,(maxCCD+1)/2); 
    }
    if (what.Contains("scope")) { c_scope= new TCanvas("c_scope","c_scope",400,800); c_scope->Divide(4,maxWF); }

    d= new DmtpcDataset;
    d->openRootFile(fname);
    for (int iccd=0; iccd<maxCCD; iccd++) {
        bias[iccd] = d->getBiasFrame(iccd+1);
	cout << "bias["<<iccd<<"] = " << bias[iccd] << endl;
    }

    // look for camera layout specification
    loadCameraLayoutMap();
}

DmtpcEventViewer::DmtpcEventViewer() {
    // Empty default constructor
    ;
}


DmtpcEventViewer::~DmtpcEventViewer() {
    // Empty desctructor
    ;
}


void DmtpcEventViewer::show(int ii) {
    iev = ii;
    validateEventNumber();
    display();
}


void DmtpcEventViewer::step(int pause_ms, int evMin, int evMax) {
  
  if (evMax == -1) evMax = tree()->GetEntries()-1;
  if (evMin == -1) evMin = 0;

  for (int ii=evMin; ii<=evMax; ii++) {
    show(ii);
    c_ccd->Update();
    gSystem->Sleep(pause_ms);
  }

}

void DmtpcEventViewer::next() {
    iev++;
    validateEventNumber();
    display();
}


void DmtpcEventViewer::prev() {
    iev--;
    validateEventNumber();
    display();
}

void DmtpcEventViewer::validateEventNumber() {

    if (iev>tree()->GetEntries()-1) {
        iev = tree()->GetEntries()-1;
        cout << "Last event" << endl;
    } else if (iev<0) {
        iev=0;
        cout << "First event" << endl;
    }
  
}

void DmtpcEventViewer::dark() {
  Int_t ipad = 0;
  
  for (int iccd=0; iccd<d->event()->ccdData()->GetEntries(); iccd++) {
    if (verbose > 1) cout << "iccd = " << iccd << endl;
    
    // determine where on the canvas to display the image
    ipad = getPadNumber(iccd);
    if (verbose > 1) cout << "ipad = " << ipad << endl;
    c_ccd->cd(ipad);
    
    // get the image from the dataset
    /* FIXME:  can you get the serial number for the dark frame? */

    if (verbose > 2) cout << "setting range: " << minRange << " to " << maxRange << endl;
    if ( (minRange==0) && (maxRange==0) ) {
      cout << "calling setFriendlyRange()" << endl;
      cout << "_nsigma = " << _nsigma << endl;
      MaxCamImageTools::setFriendlyRange(bias[iccd],_nsigma);
      cout << "bias[iccd]->GetMinimum() = " << bias[iccd]->GetMinimum() << endl;
      cout << "bias[iccd]->GetMaximum() = " << bias[iccd]->GetMaximum() << endl;
    } else {
      bias[iccd]->SetMinimum(minRange);
      bias[iccd]->SetMaximum(maxRange);          
    }
    if (verbose > 2) cout << "range set" << endl;
    
    /* FIXME:  add capability to rotate the dark frames (requires the S/N) */
    bias[iccd]->DrawCopy("colz");
    
  } 

}

void DmtpcEventViewer::displayCCD() {

  Int_t ipad = 0;
  
  for (int iccd=0; iccd<d->event()->ccdData()->GetEntries(); iccd++) {
    if (verbose > 1) cout << "iccd = " << iccd << endl;
    
    // determine where on the canvas to display the image
    ipad = getPadNumber(iccd);
    if (verbose > 1) cout << "ipad = " << ipad << endl;
    c_ccd->cd(ipad);
    
    // get the image from the dataset
    if (verbose > 2) cout << "retrieving image" << endl;
    image[iccd]=d->event()->ccdData(iccd);
    if (!image[iccd]) { cout << "image #" << iccd << " not found..." << endl; continue; }
    if (verbose > 2) cout << "image retrieved" << endl;
    //if (verbose > 1) { 
    if (verbose > 2) cout << "getting serial number for title" << endl;
    TString mytitle = TString::Format("%s evt: %d", d->event()->ccdConfig(iccd)->serialNumber.Data(), iev);
    image[iccd]->SetTitle(mytitle);
    if (verbose > 2) cout << "got serial number for title" << endl;
    //}

    if (_darkSubtract) {
      if (verbose > 2) cout << "subtracting bias" << endl;
      cout << "bias["<<iccd <<"] = " << bias[iccd] << endl;
      if (bias[iccd]) image[iccd]->Add(bias[iccd],-1);
      if (bias[iccd]) image[iccd]->SetEntries(1); // kluge
      if (verbose > 2) cout << "subtracted bias" << endl;
    }

    if (verbose > 2) cout << "setting range: " << minRange << " to " << maxRange << endl;
    if ( (minRange==0) && (maxRange==0) ) {
      cout << "calling setFriendlyRange()" << endl;
      MaxCamImageTools::setFriendlyRange(image[iccd],_nsigma);
      cout << "image[iccd]->GetMinimum() = " << image[iccd]->GetMinimum() << endl;
      cout << "image[iccd]->GetMaximum() = " << image[iccd]->GetMaximum() << endl;
      image[iccd]->SetEntries(1);
    } else {
      image[iccd]->SetMinimum(minRange);
      image[iccd]->SetMaximum(maxRange);          
    }
    if (verbose > 2) cout << "range set" << endl;
    
    if (_cameraMapExists) {
      //TH2F* tmpImg = MaxCamImageTools::rotatePerfect(image[iccd], cm->getRot(d->event()->ccdConfig(iccd)->serialNumber));
      //cout << "iccd, tmpImg = " << iccd << ", " << tmpImg << endl;
      //tmpImg->DrawCopy("colz");
      ////delete tmpImg;
      //image[iccd]->DrawCopy("colz");
      
      MaxCamImageTools::rotatePerfect(image[iccd], _cm->getRot(d->event()->ccdConfig(iccd)->serialNumber))->DrawCopy("colz");
      
    } else {
      image[iccd]->DrawCopy("colz");
    }
    
    if (what.Contains("verb")) {
      MaxCamImageTools::applyThreshold( image[iccd], 30);
      MaxCamImageTools::killLonePixels( image[iccd], 30, 4);
      cout << "Image " << iccd << " yield = " << image[iccd]->Integral() << endl;
    }
  } 
}

void DmtpcEventViewer::displayScope() {
  int nwf=d->event()->scopeData()->GetEntries()/(_nscope*2);
  int itr=0;
  for (; itr<nwf; itr++) {
    if (itr>maxWF) break;
    for (int iboard=0; iboard<_nscope; iboard++) {
      for (int ich=0; ich<2; ich++) {
	TH1F *hwf=d->event()->scopeData(itr, iboard, ich);
	hwf->SetLineColor(iboard*2+ich+1);
	hwf->SetLineWidth(3);
	c_scope->cd(itr*_nscope*2 + iboard*2+ich +1);
	hwf->Draw();
	//hwf->SetMaximum( d->event()->scopeDataInfo(i)->getVoltageMax() );
	if (what.Contains("verb")) {
	  MaxCamWaveformTools wfA(hwf); wfA.print();
	}
	
      }
    }
  }
  if (itr<nwf) cout << "(shown "<<itr << " of "<< nwf << " trigers)"<<endl;
}

void DmtpcEventViewer::display() {
  // assumes that iev is a valide event number (see validateEventNumber())
  
  cout << "Display event no. " << iev << "  with ID = " << d->event()->eventNumber() << endl;
  d->getEvent(iev);
  if (verbose > 0) cout << "after d->getEvent()" << endl;
  
  if (what.Contains("ccd")) displayCCD();

  if (what.Contains("scope")) displayScope();
    
  if (what.Contains("mc")) {
    int ntr=d->event()->mcTrack()->GetEntries();
    if (ntr>10) ntr=10;
    for (int i=0; i<ntr; i++) d->event()->mcTrack(i)->print();
    
    
    cout << "(shown " << ntr << " of " << d->event()->mcTrack()->GetEntries() << " tracks)" << endl;
  }
  
  cout << "Display event no. " << iev << "  with ID = " << d->event()->eventNumber() << endl;
  
}


TTree*
DmtpcEventViewer::tree() { return d->tree(); }

void
DmtpcEventViewer::getEvent(int iev) { d->getEvent(iev); }

DmtpcDataset*
DmtpcEventViewer::data() { return d; }

void
DmtpcEventViewer::loadCameraLayoutMap() {
  // look for a camera layout map file in the current directory
  // if none exists, do nothing else and return 
  // otherwise, parse the file and set the _cameraMapExists flag to true

  FileStat_t junk;
  if (gSystem->GetPathInfo(_cameraMapName.Data(), junk)) {  // file does not exist
    cout << "No map file found named: [" << _cameraMapName << "]" << endl;
    return;
  } else { // file exists, so read and parse the file
    cout << "Camera layout map file found: " << _cameraMapName << endl;
    _cm = new DmtpcCameraMap();
    _cm->setVerbose(2);
    _cm->loadMap(_cameraMapName);
    _cameraMapExists = true;
  }
  return;
}

Int_t 
DmtpcEventViewer::getPadNumber(Int_t iccd) {
  cout << "getPadNumber:  maxCCD = " << maxCCD << endl;
  if (maxCCD == 1) return 0;
  if (!_cameraMapExists) return iccd+1;

  cout << "d->event()->ccdConfig(iccd)->serialNumber = " << d->event()->ccdConfig(iccd)->serialNumber << endl;
  
  // otherwise, get the pad number from the Camera Layout Map
  return _cm->getPad(d->event()->ccdConfig(iccd)->serialNumber);
}

TString DmtpcEventViewer::sn(int iccd) {
  return d->event()->ccdConfig(iccd)->serialNumber;
}

void DmtpcEventViewer::drawAnodeBoundary() {

  TArc arc(0,0,1024);
  arc.SetPhimin(0);
  arc.SetPhimax(90);
  arc.SetFillStyle(0);
  arc.SetLineColor(kRed);
  arc.SetLineWidth(2);
  
  int ipad;
  for (int iccd=0; iccd<maxCCD; iccd++) {
    cout << "iccd = " << iccd << endl;
    ipad = getPadNumber(iccd);
    
    cout << "  ipad = " << ipad << endl;
    if (ipad == 1) {
      cout << "   found 1" << endl;
      arc.SetTheta(90);
      arc.SetX1(1024);
      arc.SetY1(0);
    } else if (ipad==2) {
      cout << "   found 2" << endl;
      arc.SetTheta(0);
      arc.SetX1(0);
      arc.SetY1(0);
    } else if (ipad==3) {
      cout << "   found 3" << endl;
      arc.SetTheta(180);
      arc.SetX1(1024);
      arc.SetY1(1024);
    } else if (ipad==4) {
      cout << "   found 4" << endl;
      arc.SetTheta(270);
      arc.SetX1(0);
      arc.SetY1(1024);
    }

    cout << "c_ccd->cd("<<ipad<<") " << endl;
    c_ccd->cd(ipad);
    cout << "arc.DrawClone()" << endl;
    arc.DrawClone("only");
  }


}
