#include "MaxCamRead.hh"
#include "MaxCamConfig.hh"
#include "MaxCamChannel.hh"
#include "MaxCamImageTools.hh"
#include "MaxCamMC.hh"
#include "MaxCamTrack.hh"
#include "MaxCamSegment.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
using std::cout;
using std::endl;
using std::cerr;
using std::vector;

#include "math.h"

#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TPolyMarker.h"
#include "TF1.h"
#include "TPad.h"
#include "TSystem.h"
#include "TCutG.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include "fitsio.h"

#include<map>
using std::map;
using std::pair;


MaxCamMC *MaxCamRead::srim;

double calcStopping(double *var, double *par) {
        double x  = var[0]; // x coordinate along track
        double x0 = par[0]; // left endpoint position
        double x1 = par[1]; // right endpoint position
        double calib = par[6]; // conversion pixel -> mm
        double ret   = par[3]; // background level;
        double sign  = par[5];
        if (x>x0 && x<x1) {
                double dx = sign<0? x-x0 : x1-x; 
                ret +=  par[2]*MaxCamRead::srim->getStoppingVsRange()->Eval( dx*calib );
        }
        return ret;
}

double calcStoppingConv(double *var, double *par) {
        double sigma=par[4]; // mm
        double x0=var[0];
        double xmin  = var[0]-3*sigma; // x coordinate along track
        double xmax  = var[0]+3*sigma; // x coordinate along track
        double dx=(xmax-xmin)/100;
        double ret=0;
        for (double x=xmin; x<xmax; x+=dx) {       
                var[0]=x;
                ret += TMath::Gaus(x,x0,sigma,true) * calcStopping(var, par);
        }
        var[0]=x0;
        return ret*dx;
}


ClassImp(MaxCamRead)

//____________________
// A class for reading files that contain CCD images,
// camera configuration, time stamps... that saved in a 
// ROOT tree.
// All information can be accessed event-by-event or 
// projected for all events using standard ROOT classes.
//
// Example:
//  MaxCamRead *analysis = new MaxCamRead("root://mitbbr11//data/ccdrun_00001.root");
//  analysis->event(12); // access a specific event
//  analysis->ccdImage()->Draw("colz"); // print image for given event
//  analysis->timeStamp()->Get(); // time-stamp for given event
//
//  analysis->tree()->Draw("ccdImage.Integral():timeStamp.Get()"); // plots total sum of pixels vs time
//
    MaxCamRead::MaxCamRead(TString inputFile, TString outfile) :
        MaxCamDataset( (const char *)inputFile ), 
        _biasFrame(0), _x0(-1),_x1(-1),_y0(-1),_y1(-1), _verbose(0), _ndischarge(0) {
        // A simple constructor with file name as arrument.
        // The file is opened in read-mode.
        //
        
        _fun = new TF1("_fun","[0]+[1]*exp(-0.5*(x-[2])**2/[3]**2)"); 
        
        if (outfile!="") _outputFile=new TFile(outfile,"recreate");
        nt = new TTree("res","");
        nt->Branch("setup", &Wires,     
                   "time/I:wirehv/F:meshhv:pressure:setpress:expotime/I:event:run:fitno:emc/F");
        
        _segment = new MaxCamSegment;
        nt->Branch("segment", "MaxCamSegment", &_segment, 32000, 0);
        
        wireMinROI=-1;
        wireMaxROI=-1;
        
        // zero time
        getEvent(0);
        _T0=*timeStamp();        
        
        // so, so ugly - one day I'll pay for this.
        //srim=new MaxCamMC;
        //srim->fillSrimTable("SRIM_He_in_CF4_100Torr");        
        //_fStopping = new TF1("_fStopping", calcStoppingConv, 0, 768, 7);
        //_fStopping->SetParLimits(5, -1, -1);  // negative direction  

        _recoil = 0;
        if (chain()->GetBranchStatus("recoil")) {
            _recoil = new TLorentzVector;
            chain()->SetBranchAddress("recoil", &_recoil);
            cout << "Found Recoil Lorentz Vector for recoil" << endl;
             nt->Branch("mcrecoil", "TLorentzVector", &_recoil, 32000, 0);
        }
        _projectile=0;
        if (chain()->GetBranchStatus("projectile")) {
            _projectile = new TLorentzVector;
            chain()->SetBranchAddress("projectile", &_projectile);
            cout << "Found Recoil Lorentz Vector for projectile" << endl;
             nt->Branch("mcprojectile", "TLorentzVector", &_projectile, 32000, 0);
        }
        
}


void
MaxCamRead::addRun(TString fname) {
    // Add run file to the current collection of events.
    chain()->Add(fname);
}

int
MaxCamRead::getRunNumber() {
        TString runname(chain()->GetCurrentFile()->GetName());
        int first=runname.First('_')+1;
        int runno=atoi( &runname[first]);
        Wires.run = runno;
        return runno;
}

void
MaxCamRead::setEventInfo(int iev) {
        Wires.event=iev;
        // some basic event info
        Wires.time   =
            (timeStamp()->GetDay()  - _T0.GetDay())*86400 +
            (timeStamp()->GetHour() - _T0.GetHour())*3600 +
            (timeStamp()->GetSecond() - _T0.GetSecond());
        Wires.wirehv = wire()->setValue;
        Wires.meshhv = mesh()->currentValue;
        Wires.pressure = pressure()->currentValue;  
        Wires.setpress = pressure()->setValue;  
        Wires.expotime = ccdConfig()->exposureTime;
        Wires.emc = _recoil ? _recoil->P() : 0;
}


TH1F* 
MaxCamRead::makeYieldHisto(TH2F *img, float min, float max) {
  // Make 1D histogram of pixel yields from current image.
  //

    TH2F* tmpImage = img ? img : ccdImage();

  if (min<1 || min>max) min = tmpImage->GetMinimum();
  if (max<1 || min>max) max = tmpImage->GetMaximum();
  int nbin = int(max-min);

  TH1F* hY=new TH1F("hY","", nbin, min,max);

  if (_x0<0 || _x1<0 || _y0<0 || _y1<0) {
    _x0=1; _x1=tmpImage->GetXaxis()->GetNbins();
    _y0=1; _y1=tmpImage->GetYaxis()->GetNbins();
  }

  //cout <<"x0,x1="<<_x0<<","<<_x1<<"   y0,y1="<<_y0<<","<<_y1<<endl;

  for (int i=_x0; i<=_x1; i++) {
    for (int j=_y0; j<=_y1; j++) {
      //cout << tmpImage->GetBinContent(i,j)<< " ";
      hY->Fill( tmpImage->GetBinContent(i,j) );
    }
  }
  return hY;
}

void 
MaxCamRead::setROI( int x0, int x1, int y0, int y1) {
  // Set ROI for display.
  //

  _x0=x0; _x1=x1; _y0=y0; _y1=y1;
  ccdImage()->GetXaxis()->SetRange(_x0,_x1);
  ccdImage()->GetYaxis()->SetRange(_y0,_y1);
}

void
MaxCamRead::getEvent(int i) {
  // Access event tree for a specific event.
  //
  _imageTree->GetEvent(i); 
  
  if (_biasFrame) ccdImage()->Add(_biasFrame,-1);

  for (vector<int>::iterator ih=_hotPixels.begin(); ih!=_hotPixels.end(); ih++) {
      ccdImage()->SetBinContent(*ih,0);
  }
}

void
MaxCamRead::readBiasFrame() {
  // Check if bias frame was saved in the ROOT file and correct
  // images substracting the bias frame. If the bias frame was not
  // found no correction is applied.
  
  deleteBiasFrame();

  _biasFrame = (TH2F*)_file->Get("bias");
}

void
MaxCamRead::deleteBiasFrame() {
  // Remove bias frame from analysis - from now on event images
  // do not have automatic ADC bias frams substraction.

  if (_biasFrame) {
    delete _biasFrame;
    _biasFrame=0;
  }  
}




void
MaxCamRead::findHotPixels(TString hotFile) {
    if (_verbose) cout << endl << "Find hot pixels: " << hotFile << endl;;
  deleteHotPixels();

  ifstream hf((const char*)hotFile);
  int ipix;
  while (!hf.eof()) {
    hf >> ipix;
    _hotPixels.push_back( ipix );
    //cout << ipix << endl;
  }
}



void
MaxCamRead::findHotPixels(int nevents, float th, float rate) {

    deleteHotPixels();
    
    if (nevents>tree()->GetEntries()) nevents=tree()->GetEntries();
    
    // find pixels over 10 sigma in nevents
    map<int,int> hotMap;
    for (int i=0; i<nevents; i++) {
        getEvent(i);
        float mean = ccdImage()->GetMean();
        float rms  = ccdImage()->GetRMS();
        int nx = ccdImage()->GetNbinsX();
        int ny = ccdImage()->GetNbinsY();
        for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
                if (ccdImage()->GetBinContent(i,j)-mean<rms*th) continue;
                map<int,int>::iterator ih = hotMap.find( ccdImage()->GetBin(i,j) );
                if ( ih==hotMap.end()) hotMap.insert( pair<int,int>( ccdImage()->GetBin(i,j),1) );
                else ih->second++;
            }
        }
        
    }
    
    // find those repeating 
    for (map<int,int>::iterator ih = hotMap.begin(); ih!=hotMap.end(); ih++) {
        if (ih->second<nevents*rate)  continue;
        _hotPixels.push_back( ih->first );
        if (_verbose) cout << "hot pixel: " << ih->first << "  rate=" << float(ih->second)/nevents<<endl;
    }
    
}

void
MaxCamRead::deleteHotPixels() {
  _hotPixels.clear();
}


void
MaxCamRead::findWires(float width, float thr, int nev, TString opt, float maxAvgSum) {
  // Makes a search for wires by accumulating many events.
  // It assumes wires are along y or x direction so a projection of
  // accumulated images contains peaks of light collected
  // by different wires. It filles wireList vector with
  // coordinates of wire positions. 'width' and 'threshold' are parameters for
  // for TH1F:ShowPeaks, 'nev' is the number of events used in the search (<0 means all), 
  // 'maxAvgSum' is the maximum sum of all pixels (<0 means no cut),
  // 'opt' gives the orientation of wires (x or y)
  //


  // add frames
  int n= nev>0 ? nev : tree()->GetEntries();
  getEvent(0);
  _xproj = opt.Contains("x") ? 
    (TH1D*)ccdImage()->ProjectionX()->Clone("_xproj") :
    (TH1D*)ccdImage()->ProjectionY()->Clone("_xproj");
  
  for (int ii=1; ii<n; ii++) {
    getEvent(ii);
     if (maxAvgSum>0 && ccdImage()->Integral()>maxAvgSum) continue;
    _xproj->Add( opt.Contains("x") ? ccdImage()->ProjectionX() : ccdImage()->ProjectionY()  );
  }
  // search for peaks
  wireList.clear();  
  _xproj->ShowPeaks(width, "", thr);
  TPolyMarker *wireX = (TPolyMarker*)_xproj->GetFunction("TPolyMarker");
  for (int ii=0; ii<wireX->GetN(); ii++) wireList.push_back( wireX->GetX()[ii] );
  std::sort(wireList.begin(), wireList.end());
  createWireBinList();

}


TH1D* 
MaxCamRead::wireYield(int iwire, int dwire, int type, TString opt) { 
    // Evaluate light intensity at given wire (assumes 1-alpha track per event)
    // and fill Wire structure that is saved into an ntuple.
    // Options are:
    //    p = plot each fit (pause for 1sec)
    //    x = wire defined by x coordinate
    //    y = wire defined by y coordinate
    //    
    // Type parameter:
    //    0 = fit signal along wire ('dwire' is the width of the strip)
    //    1 = fit signal perp. to wire (e.g. to get integral of all hits on the wire)
    //

  TH1D *hy=0;

  Wire[iwire].y = Wire[iwire].yerr = Wire[iwire].ampl = 0; 
  Wire[iwire].mean = Wire[iwire].width = Wire[iwire].chi2 = 0;
  Wire[iwire].ndof = 0;
  Wire[iwire].id = iwire;

  if (opt.Contains("0")) return hy; // option to skip fit

    
  // transverse ROI
  int xwire=wireBinList[iwire];
  int wireMinPixel=xwire-dwire;
  int wireMaxPixel=xwire+dwire;


  // find intensity
  switch (type) {
  case 0: 
      if (opt.Contains("x")) hy=ccdImage()->ProjectionY("_py", wireMinPixel, wireMaxPixel); 
      else                   hy=ccdImage()->ProjectionX("_px", wireMinPixel, wireMaxPixel);
      break;
  case 1: 
      if (wireMinROI>=wireMaxROI || wireMinPixel>=wireMaxPixel) return hy;
      if (opt.Contains("x")) hy=ccdImage()->ProjectionY("_px", wireMinROI, wireMaxROI); 
      else                   hy=ccdImage()->ProjectionX("_py", wireMinROI, wireMaxROI); 
      hy->GetXaxis()->SetRange(wireMinPixel,wireMaxPixel);
      break;
  }


  
  // find bin density
  double brho = hy->GetNbinsX()/(hy->GetXaxis()->GetXmax()-hy->GetXaxis()->GetXmin());

  double bkgr = hy->GetMinimum();
  double ampl = hy->GetMaximum() - hy->GetMinimum();
  double mean = hy->GetBinCenter( hy->GetMaximumBin() );
  double width= 10;
  if (type==1) {
      mean= wireList[iwire];
      width=5;
  }
  _fun->SetParameter(0, bkgr);
  _fun->SetParameter(1, ampl);
  _fun->SetParameter(2, mean);
  _fun->SetParameter(3, width);

  
  //_fun->SetParLimits(2, mean, mean);
  //_fun->SetParLimits(3, width, width);
  //_fun->SetParLimits(3, 3, 13);

  //_fun->SetParLimits(2, mean-10, mean);

  if (opt.Contains("p")) {
      hy->Fit("_fun","Q"); hy->Draw(); gPad->Update(); gSystem->Sleep(1000);
  } else hy->Fit("_fun","Q0");


  static double sqrt2pi=sqrt(2*TMath::Pi());
  float I = hy->GetFunction("_fun")->GetParameter(1)
      * fabs( hy->GetFunction("_fun")->GetParameter(3) )
      * sqrt2pi * brho;

  
  float IErr = I* sqrt( pow(hy->GetFunction("_fun")->GetParError(1)/hy->GetFunction("_fun")->GetParameter(1),2) +
			pow(hy->GetFunction("_fun")->GetParError(3)/hy->GetFunction("_fun")->GetParameter(3),2) );


  // save fit
  Wire[iwire].y     = I; 
  Wire[iwire].yerr  = IErr;
  Wire[iwire].ampl  = hy->GetFunction("_fun")->GetParameter(1); 
  Wire[iwire].mean  = hy->GetFunction("_fun")->GetParameter(2);
  Wire[iwire].width = hy->GetFunction("_fun")->GetParameter(3);
  Wire[iwire].chi2  = hy->GetFunction("_fun")->GetChisquare();
  Wire[iwire].ndof  = hy->GetFunction("_fun")->GetNDF();
  Wire[iwire].id    = iwire;

  if (_verbose) {
      cout << "I = " << I << " +- " << IErr << endl;
  }


  return hy; 
}



TH1F*
MaxCamRead::mergeWireStrips(int dwire0, int dwire1, TString opt, TString hname) {
    // Create a projection of pixels that are within range (dwire0,dwire1) 
    // from each wire. Background is not subtracted.
    
    TH2F* image = opt.Contains("c") ? ccdImage() : _imageth;
    
    bool isX = opt.Contains("x");
    
    int nbin = isX ? image->GetNbinsY() : image->GetNbinsX();
    TH1F *project = new TH1F(hname, "", nbin, 0, image->GetXaxis()->GetXmax() );
    int nperpbin = isX ? image->GetNbinsX() : image->GetNbinsY();
    
    //cout << "WireIdList size is " << wireIdList.size() << endl;
    
    for (unsigned int iw=0; iw<wireIdList.size(); iw++) {
        
        int wireMinPixel=wireBinList[ wireIdList[iw] ]+dwire0;
        int wireMaxPixel=wireBinList[ wireIdList[iw] ]+dwire1;
        if (_verbose) cout << "From: " << wireMinPixel << " - " << wireMaxPixel << endl;
        if (wireMinPixel<1 || wireMaxPixel>nperpbin) {
            cout << GetName() << ": selected region over bounds! Asserting..."<<endl;
            assert(0);
        }
        
        // fill histogram projection with threshold cut
        float y=0;
        for (int i=wireMinPixel; i<=wireMaxPixel; i++) {
            for (int j=1; j<=nbin; j++) {
                y = isX ? image->GetBinContent(i,j) : image->GetBinContent(j,i); 
                project->SetBinContent( j, project->GetBinContent(j)+y);
            }
        }
    }
    return project;
}



TH1F*
MaxCamRead::findSegmentLength(int iwire, int dwire, float threshold, int pixelGapMaximum, TString opt) {
    // Numeric approach to find segment length at a wire at low occupancies.
    // It merges all wires into one projection and
    // subtracts background (see options in makeWireProjection)
    // finally, finds the length of the segment (see findSegmentLength)
    
    
    // make image with threshold
    _imageth = (TH2F*)ccdImage()->Clone("imageth");
    MaxCamImageTools::applyThreshold(_imageth, threshold);
    MaxCamImageTools::killLonePixels(_imageth, threshold);
    
    // project wires that are in wireIdList
    wireIdList.clear();
    wireIdList.push_back(iwire);
    TH1F *project = mergeWireStrips( -dwire,    dwire, opt, "project");
    assert(project);
    int iR, iL;
    
    if (opt.Contains("s")) { // fit to stopping power function
        project->Fit("_fStopping");
        Wire[iwire].widthL=project->GetFunction("_fStopping")->GetParameter(0);
        Wire[iwire].widthR=project->GetFunction("_fStopping")->GetParameter(1);
        Wire[iwire].ycount=project->GetFunction("_fStopping")->GetParameter(2);
        Wire[iwire].shapeChi2=project->GetFunction("_fStopping")->GetChisquare(); 
        iL=project->FindBin( Wire[iwire].widthL );
        iR=project->FindBin( Wire[iwire].widthR );
        project->GetXaxis()->SetRange( iL, iR );
        Wire[iwire].skewness=project->GetSkewness(1);
    }

    else { // simple counting

        // get bkg level as mean + rms of pixel yields in the projection of background (bg)
        int dw0=-dwire*3;
        int dw1=-dwire;
        if (wireBinList[ iwire ] + dw0<0) {
            dw0=dwire;
            dw1=3*dwire;
        }
        
        TH1F *projectBkg   = mergeWireStrips( dw0, dw1, opt, "sideband");
        TH1F *projectYield = MaxCamImageTools::createYieldHisto( projectBkg );
        double bkg0 = projectYield->GetMean();
        double bkg  = bkg0 + projectYield->GetRMS();
        delete projectYield;
        delete projectBkg;
        
        MaxCamImageTools::countSegmentLength(project, bkg, pixelGapMaximum, iR, iL);
        Wire[iwire].widthL = float(iL);
        Wire[iwire].widthR = float(iR);
        Wire[iwire].skewness = MaxCamImageTools::calcSkewness( project, int(Wire[iwire].widthL), int(Wire[iwire].widthR), bkg0);
    }
    
    
    double rangeInmm = ( Wire[iwire].widthR - Wire[iwire].widthL )/(96./12.) ;
    TGraph *range = srim->getEnergyVsRangeProject( "long" );
    Wire[iwire].esrim = range->Eval( rangeInmm );
    delete range;
    
    if (_verbose) {
        cout <<"Edges = " << Wire[iwire].widthL << ",  " << Wire[iwire].widthR << endl;
        cout << "Skewness=" << Wire[iwire].skewness << endl;
        cout << "Range=" << rangeInmm 
             << "   Energy=" << Wire[iwire].esrim   << endl;
    }
    
    delete _imageth;
    return project;
}




void 
MaxCamRead::makeNtuple(int minwire, int maxwire, int width, int type, TString opt) {
   // Assuming wires are parallel with x or y axis (specify direction as input option), loop over events
   // and find hits on the wires. Each hit is assumed to be a gaussian, and only one hit per wire 
   // is allowed. A ntuple is created with hit positions, widths, integrals for each wire in range
   // 'minwire' and 'maxwire' in the wire list (see findWires for creating a list of wires). 
   // 'width' specifies the width of the wire in ccd pixels.
   // 'type'=0 makes a projection along wire, 1 makes projection perp. to wire direction.
   // If option string contaons 'p' it will display each fit fro 1sec, which is useful for
   // checking that fits make sense.
   //
   //
   // Example:
   //    MaxCamRead *ana=new MaxCamRead("testfile.root");
   //    ana->readBiasFrame();
   //
   //    ana->findWires(4,0.1,"y"); // assume wires are || with y-axis
   //    ana->makeNtuple(0,2,3,0,"y"); // make ntuple from hits on wires 0,1,2
   //

    int n=tree()->GetEntries();
    getEvent(0);
    //_time0=timeStamp()->Get();
    //_day0 =timeStamp()->GetDay();


    for (int i=0; i<n; i++) {	
    	if (_verbose) cout <<"event="<<i<<endl;
    	getEvent(i);
    	for (int ii=minwire; ii<=maxwire; ii++) { 
	    wireYield(ii, width, type, opt); 
	    nt->Fill(); 
	}
    }  
}



int
MaxCamRead::convertIntoFits() {
  // Convert all image histograms in the current ROOT file into
  // FITS images. Each image will be saved into a separate file
  // that will contain the base name of the ROOT file and the image
  // number. See also MaxCamImageTools.
  //
  // Example:
  //   MaxCamRead *analysis = new MaxCamRead("root:://mitbbr11//data/ccdrun_00001.root");
  //   analysis->convertIntoFits();
  //

  int n=tree()->GetEntries();
  TString fNameBase=_file->GetName();
  fNameBase.ReplaceAll(".root","_");
  fNameBase.ReplaceAll("root://","");
  fNameBase.ReplaceAll("/","_");
  
  int status = 0;
  for (int ii=0; ii<n; ii++) {
    tree()->GetEvent(ii);

    //long naxes[2] = { ccdConfig()->row_width, ccdConfig()->img_rows};
    //fitsfile *fp;

    TString fileName = fNameBase + Long_t(ii) + ".fits";
    status=MaxCamImageTools::convertIntoFits( ccdImage(), fileName);
    if (status) break;
  }

  
  TString biasName = fNameBase + "_bias.fits";
  readBiasFrame();
  status=MaxCamImageTools::convertIntoFits( _biasFrame, biasName);
  
  
  return status;
}



TH2F*
MaxCamRead::accumulatePressure(float press, float avgmin, float avgmax, int fst, int nev) {
    // Select all images with given pressure and add them up
    // into a single image.
    
    TH2F *hall=0;
    int n0 = fst>0 ? fst : 0; 
    int n = nev>0 ? nev : tree()->GetEntries(); 
    for (int i=n0; i<n; i++) {
	getEvent(i);

        float avg=ccdImage()->Integral()/(ccdImage()->GetNbinsX()*ccdImage()->GetNbinsY());
        if( avg<avgmin || avg>avgmax) continue;
        
	if (press>0 && fabs(pressure()->currentValue-press)>1e-3) continue;
        
	if (!hall) hall=(TH2F*)ccdImage()->Clone("hall");
	else hall->Add( ccdImage() );
    }
    return hall;
}



void 
MaxCamRead::select(TH2F* image, TString opt) {
    // Select region defined by '*cutSG', a TCutG object.
    
    int nx=image->GetNbinsX();
    int ny=image->GetNbinsY();

    for (int i=1; i<=nx; i++) {
	float xbin = image->GetXaxis()->GetBinCenter(i);
	for (int j=1; j<=ny; j++) {
	    float ybin = image->GetYaxis()->GetBinCenter(j);
	    bool pass = opt.Contains("sg") ? _cutSG->IsInside(xbin,ybin) : _cutBG->IsInside(xbin,ybin);
	    if (pass) continue;
	    image->SetBinContent(i,j,0);
	}
    }
}


void
MaxCamRead::setRectangularCut(float x0, float y0, float x1, float y1, float d, TString opt) {
    // Set a rectangular cut with one axis defined as (x0,y0)-(x1,y1) and the other
    // perpendicular axis with depth d.

  _xaxis = new TVector2(x1-x0, y1-y0);
  *_xaxis = _xaxis->Unit();
  _yaxis = new TVector2( -_xaxis->Y(), _xaxis->X() );

  float x[]={ x0, x1, x1+d*_yaxis->X(), x0+d*_yaxis->X(), x0};	
  float y[]={ y0, y1, y1+d*_yaxis->Y(), y0+d*_yaxis->Y(), y0};
  if (opt.Contains("sg")) { _cutSG = new TCutG("cutSG",5,x,y); }
  else  { _cutBG = new TCutG("cutBG",5,x,y); _cutBG->SetLineStyle(2); _cutBG->SetMarkerStyle(24); }
}



TH2F* 
MaxCamRead::rotate(TH2F* image, int nrot) {
    // Rotate the image using axes defined by a rectangular cut region (see setRectangularCut 
    // for defining the cut).

    double x0, y0;
    _cutSG->GetPoint(0,x0,y0);

    return (TH2F*) MaxCamImageTools::rotate(image, nrot, _xaxis, _yaxis, x0, y0);
}



  
void 
MaxCamRead::addHotPixel(int binx, int biny) { 
    // Add a chanel to hot pixel list

    if (!ccdImage()) {
        cout << GetName() << ": image not found - cannot add hot pixel" << endl;
        return;
    }
    int bin=ccdImage()->GetBin(binx, biny);
    cout << GetName() << ": adding channel " << bin << endl; 
    _hotPixels.push_back( bin);
}


void
MaxCamRead::createWireBinList() {
    wireBinList.clear();
    for (unsigned int ii=0; ii<wireList.size(); ii++) {
        wireBinList.push_back( _xproj->FindBin( wireList[ii] ) );
        cout << "wire " << ii << "= " << wireList[ii] << "   " << wireBinList[ii] << endl;
    }
}


void
MaxCamRead::addWire(double xwire) {
    // Add wire to the wire list.
    
    wireList.push_back( xwire );
    std::sort(wireList.begin(), wireList.end());
    createWireBinList();
}


void
MaxCamRead::removeWire(int iwire) {

    wireList.erase( (vector<double>::iterator)&wireList[iwire] );
    createWireBinList();
}



bool
MaxCamRead::isDischargeEvent(float threshold) {
    bool ret=
        ccdImage()->Integral()
        > threshold*ccdImage()->GetNbinsX()*ccdImage()->GetNbinsY();
    if (ret) _ndischarge++;
    return ret;
}


bool
MaxCamRead::hasTracks() {
    bool ret =false;
    TH1F *hY = makeYieldHisto();
    MaxCamTrack* trfit = new MaxCamTrack( ccdImage(), true );
    float threshold = hY->GetMean() + hY->GetRMS()*3;
    trfit->setThreshold( threshold );
    trfit->makeTracks();
    if (trfit->nTracks()) { ret=true; }
    delete hY;
    delete trfit;
    if (ret) return true;
    
    /* trfit = new MaxCamTrack( ccdImage(), false );
    trfit->setThreshold( threshold );
    trfit->makeTracks();
    if (trfit->nTracks()) { ret = true; }
    delete hY;
    delete trfit;
    */
    return ret;
}


bool
MaxCamRead::hasSegments(int &imax, int &jmax, float &threshold) {
    bool ret = false;
    imax=jmax=0;
    
    TH1F *hY = makeYieldHisto();
    threshold = hY->GetMean() + hY->GetRMS()*3;
    delete hY;

    
    MaxCamImageTools::applyThreshold(ccdImage(), threshold);
    MaxCamImageTools::killLonePixels(ccdImage(), threshold);
    vector<int> bins=MaxCamImageTools::killSecondaryClusters(ccdImage(), threshold);
    cout << "BINS SIZE="<< bins.size()<<endl;
    MaxCamImageTools::countSegmentLength2D(ccdImage(), threshold, imax, jmax);
    if (imax && jmax) {
        ret=true;

        // check if segment touching boundary
        //vector<int> bins=MaxCamImageTools::findNeighbors(ccdImage(), threshold, imax);
        for (unsigned int i=0; i<bins.size(); i++) {
            int xi = bins[i]%(ccdImage()->GetNbinsX()+2);
            int yi = bins[i]/(ccdImage()->GetNbinsX()+2);
            if (xi<3 || yi<3 || xi>ccdImage()->GetNbinsX()-2 || yi>ccdImage()->GetNbinsY()-2) {
                segment()->setBorderline(true);
                break;
            }
        }

        // number of pixels in segment
        segment()->setNPixels( bins.size() );

        // find segment length = maximumseparation between pixels
        segment()->setLength( MaxCamImageTools::calcPixelDistance(ccdImage(), imax, jmax), 0 );

        // find correlation coef. (->change to find moments)
        segment()->setCorrelation( MaxCamImageTools::calcPixelCorrelation(ccdImage(), imax, jmax) );

        // principal moments of inertia
        float Ix, Iy;
        MaxCamImageTools::principalAxes( ccdImage(), threshold, Ix, Iy);
        segment()->setIxy( Ix, Iy);

        // cos recoil 2D
        segment()->setCosRecoil (MaxCamImageTools::cosRecoil2D( ccdImage(), imax, jmax ) );
        
        
    }
    bins.clear();    
    
    return ret;
}



