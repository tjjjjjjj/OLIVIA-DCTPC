#include "MaxCamTrack.hh"

#include "MaxCamImageTools.hh"

#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TVector2.h"

#include <math.h>
#include <iostream>
using std::cout;
using std::endl;

MaxCamTrack::MaxCamTrack() {;} 

MaxCamTrack::MaxCamTrack(TH2F *img, bool rotate90deg) : 
  TObject(),
  _th(100),
  _wireCoord("y"),
  _pihalf(TMath::Pi()*0.5),
  _yroi(400),
  _phiroi(5)
{
  
  if (rotate90deg) {
      //cout << "rotating image" << endl;
      _image = new TH2F("_image","",
                        img->GetNbinsY(), img->GetYaxis()->GetXmin(), img->GetYaxis()->GetXmax(),
                        img->GetNbinsX(), img->GetXaxis()->GetXmin(), img->GetXaxis()->GetXmax() );
      for (int i=1; i<=img->GetNbinsX(); i++)
          for (int j=1; j<=img->GetNbinsY(); j++) 
              _image->SetBinContent(j,i, img->GetBinContent(i,j) );
  }
  else _image = img;
  
  _nx = _image->GetNbinsX();
  _ny = _image->GetNbinsY();
  _npixel = _image->GetNbinsX() * _image->GetNbinsY();
  _xratio = _image->GetXaxis()->GetXmax()/_nx;
  _yratio = _image->GetYaxis()->GetXmax()/_ny;
  _imageth=0;
  _hphi=0;
  _hphibins=100;
  _hphimin=-_pihalf;
  _hphimax= _pihalf;
  _imagephi=0;
  _ftrack=0;
  //_neighbors=-1;
}


MaxCamTrack::MaxCamTrack(const MaxCamTrack &other) {
    // Not yet done
}


MaxCamTrack::~MaxCamTrack() {
  clearMemory();
}

void 
MaxCamTrack::checkThresholdAndNeighbors() {

  _imageth = (TH2F*)_image->Clone("imageth");

  // first apply threshold
  MaxCamImageTools::applyThreshold( _imageth, _th);

  // kill lone pixels
  MaxCamImageTools::killLonePixels(_imageth,_th);

  // kill pixels not on wires
  MaxCamImageTools::killUnusedPixels(_imageth, _wireBinList, 3, _wireCoord);
}


bool
MaxCamTrack::isNew(int ix, int iy) {
    int bin=_imageth->GetBin(ix,iy);
    int npl=_peakList.size();
    for (int i=0; i<npl; i++) {
        if (bin==_peakList[i]) return false;
    }
    return true;
}


void
MaxCamTrack::findMaxBin() {
  _i0=-1;
  _j0=-1;

  float max=_th;
  for (int i=1; i<=_nx; i++) {
    for (int j=1; j<=_ny; j++) {
      int bin = _imageth->GetBin(i,j);
      float yield = _imageth->GetBinContent(bin);

      if (yield>max && MaxCamImageTools::hasNeighbor(_imageth,i,j,_th) && isNew(i,j) ) {
	max = yield;
	_i0=i;
	_j0=j;
      }
    }
  }
  _peakList.push_back( _imageth->GetBin(_i0, _j0) );
  //cout << _i0 << "  " << _j0 << endl;
}



void
MaxCamTrack::findSlope(int ic, int jc) {

  if (ic>0 && jc>0) {  _i0=ic; _j0=jc; }
  if (_i0<1 || _j0<1) {
    if (_hphi) delete _hphi;
    _hphi=0;
    return;
  }

  //if (!_hphi) _hphi = new TH1F("hphi","", 100,-_pihalf, _pihalf);
  if (!_hphi) _hphi = new TH1F("hphi","", _hphibins, _hphimin, _hphimax);
  _hphi->Reset();

  //if (!_imagephi) _imagephi = (TH2F*)_image->Clone("imagephi");
  //_imagephi->Reset();

  for (int i=1; i<=_nx; i++) {
    for (int j=1; j<=_ny; j++) {
      if (i==_i0 && j==_j0) continue;
      if (_imageth->GetBinContent(i,j)<_th) continue;    
      float phi =  atan2(j-_j0,i-_i0);
      while (phi>_hphimax) phi-=TMath::Pi();
      while (phi<_hphimin) phi+=TMath::Pi();
      _hphi->Fill( phi, _imageth->GetBinContent(i,j)  );
      //_imagephi->SetBinContent(i,j,phi);
    }
  }

  //cout <<"_phi="<<_hphi<<endl;
}



void 
MaxCamTrack::fitTrack() {

  if (!_hphi) {
    if (_ftrack) delete _ftrack;
    _ftrack=0;
    return;
  }

  // find slope
  int maxphibin = _hphi->GetMaximumBin();
  //_hphi->GetXaxis()->SetRange( maxphibin-_phiroi, maxphibin+_phiroi); 
  float slope = tan((_hphi->GetBinCenter( maxphibin ) ) *_yratio/_xratio);

  // create formula 
  float xcoord = (_i0-0.5)*_xratio;
  float ycoord = (_j0-0.5)*_yratio;

  if (_ftrack) delete _ftrack;

  _ftrack = new TF1("ftrack","[0]+[1]*(x-[2])", 
		    _image->GetXaxis()->GetXmin(), 
		    _image->GetXaxis()->GetXmax() );

  _ftrack->SetLineStyle(3);

  _ftrack->SetParameters(ycoord, slope, xcoord);
  _ftrack->SetParLimits(2, xcoord, xcoord);

  vector<float> x, y, xerr, yerr;

  float feval;

  for (int i=1; i<=_nx; i++) {
    for (int j=1; j<=_ny; j++) {

      if ( _imageth->GetBinContent(i,j)<_th ) continue;
 
      feval = _ftrack->Eval( (i-0.5)*_xratio );
      if (fabs(feval-(j-0.5)*_yratio)>_yroi) continue;

      x.push_back( i*_xratio );
      y.push_back( j*_yratio );
      xerr.push_back( 0 );
      yerr.push_back( 1./sqrt( _imageth->GetBinContent(i,j) ) );
      
    }
  }  
  
  TGraphErrors *gr = new TGraphErrors((int)x.size(),&x[0],&y[0],&xerr[0],&yerr[0]);
  gr->Fit("ftrack","q0");
  _ftrack->SetParameter(0, gr->GetFunction("ftrack")->GetParameter(0) );
  _ftrack->SetParameter(1, gr->GetFunction("ftrack")->GetParameter(1) );
  delete gr;
}


TVector2
MaxCamTrack::distanceFromPoint(TF1 *t, double x, double y) {
  // Distance of point from the track.

  float a=t->GetParameter(1);
  float b=t->Eval(0);
  
  float newa=-1./a;
  float newb=y-x*newa;
  
  float xcross=(newb-b)/(a-newa);
  float ycross=a*xcross+b;


  return TVector2( x-xcross, y-ycross);
}



void 
MaxCamTrack::cleanUpPixels() {
  // Clean pixels around newly found track.
  //

  double delta=0;
  for (int i=1; i<=_nx; i++) {
    for (int j=1; j<=_ny; j++) {
      if ( _imageth->GetBinContent(i,j)<_th ) continue;  // dont spend time on unused pixels
      delta = distanceFromPoint( _ftrack, _image->GetXaxis()->GetBinCenter(i), _image->GetYaxis()->GetBinCenter(j) ).Mod();
      if (delta<_yroi) _imageth->SetBinContent(i,j,0);
    }
  }
  
}



void 
MaxCamTrack::makeTracks(TCanvas *dbgCanvas) {

    _peakList.clear();
    
    checkThresholdAndNeighbors();

    int _maxTrackSearchLoop=10000;
    while(_maxTrackSearchLoop--) {
        
        findMaxBin();
        if (_i0==-1 || _j0==-1) break;
       
        findSlope();        

        fitTrack();
        
        plotDebug(dbgCanvas); //return;;
        
        if (!getTrack()) continue;
        
        trackList.push_back( (TF1*)getTrack()->Clone() );

        
        cleanUpPixels();
    } 

    delete _imageth;
    _imageth=0;
    if (TString(_image->GetName())=="_image") {
        delete _image;
        _image=0;
    }

}


void
MaxCamTrack::plotDebug(TCanvas *cdbg) {
    if (!cdbg) return;
    cdbg->cd(1); _imageth->DrawCopy("colz");
    cdbg->cd(2); _hphi->DrawCopy();
    cdbg->cd(3); _image->DrawCopy("colz");
    for (unsigned int i=0; i<trackList.size(); i++) {
        trackList[i]->Draw("same");
    }
}

TH1F* 
MaxCamTrack::makeResiduals(int itr, TString opt, float width) {

  if (!trackList[itr]) {
    cout << "MaxCamTrack: track " << itr << " does not exist" << endl;
    return 0;
  }
  if (!_image) {
    cout << "MaxCamTrack: Image not set " << endl;
    return 0;
  }

  if (width<0) width = _yroi*3;
  TH1F *hres=new TH1F("hres","",2*int(width), -width, width);


  bool aboveThreshold = opt.Contains("th");
    
  // define sign of residuals
  TVector2 ortho = distanceFromPoint( getTrack(itr), 
				      _image->GetXaxis()->GetBinCenter(_nx/2), 
				      _image->GetYaxis()->GetBinCenter(_ny/2) ).Unit();  
  
  TVector2 delta;
  for (int i=1; i<=_nx; i++) {
    for (int j=1; j<=_ny; j++) {
      if (aboveThreshold && _image->GetBinContent(i,j)<_th) continue;

      delta = distanceFromPoint( getTrack(itr), 
				 _image->GetXaxis()->GetBinCenter(i), 
				 _image->GetYaxis()->GetBinCenter(j) );
      /*cout << _image->GetXaxis()->GetBinCenter(i) << "  "
	   << _image->GetYaxis()->GetBinCenter(j) << "  "
	   << delta*ortho 
	   << endl;*/
      
      hres->Fill( delta * ortho,  _image->GetBinContent(i,j) );
      
    }
  }

  return hres;
}


void
MaxCamTrack::clearMemory() {
  // Clear images used in the reconstruction from memory.

  if (_hphi) { delete _hphi; _hphi=0; }

  if (_image) { delete _image; _image=0; }

  if (_imageth) { delete _imageth; _imageth=0; }

  if (_imagephi) { delete _imagephi; _imagephi=0; }  

  for (unsigned int i=0; i<trackList.size(); i++) delete trackList[i];
  trackList.erase( trackList.begin(), trackList.end());

}
