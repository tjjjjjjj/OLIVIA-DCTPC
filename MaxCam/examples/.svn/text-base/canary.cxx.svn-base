#include "../DmtpcDataset.hh"
#include "../MaxCamImageTools.hh"
#include "../DmtpcEvent.hh"
#include "../MaxCamSRIM.hh"
#include "../MaxCamTrack.hh"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TLegend.h"
//#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"
#include "TSystem.h"

#include <vector>
using std::vector;

#include <iostream>
using std::cout;
using std::endl;
using std::flush;

DmtpcDataset *d=0;
TTree *nt=0, *ft=0;

float sqrt2pi=sqrt(TMath::TwoPi());

int iCCD=0;





//////////////////////////////////////////////////////////////////////////
//
//	 MAKE BIAS FRAME
//
//////////////////////////////////////////////////////////////////////////



TH2F *bias=0;
void makeBiasFrame(TString fname) {
    if (bias) { delete bias; bias=0; }
    d= new DmtpcDataset;
    d->openRootFile(fname);

    int nev=d->tree()->GetEntries();
    for (int i=0; i<nev; i++) {
        d->getEvent(i);
        if (!bias) bias=(TH2F*)d->event()->ccdData(iCCD)->Clone("bias");
        else bias->Add( d->event()->ccdData(iCCD) );
    }
    bias->Scale(1./nev);
    delete d; d=0;
}


//////////////////////////////////////////////////////////////////////////
//
//	NOISE STUDY
//
//////////////////////////////////////////////////////////////////////////

vector<float> ccdReso, expot;
void ccdResolution(TString fname, TString opt="") {
	if (!bias) makeBiasFrame(fname);
	DmtpcDataset*  dd = new DmtpcDataset;
	dd->openRootFile(fname);
	dd->getEvent(1);
	TH2F *image = (TH2F*)dd->event()->ccdData(iCCD)->Clone("image");
	image->Add(bias,-1);
	TH1F * hist=MaxCamImageTools::createYieldHisto( image, -2000, 3000 );
	hist->Fit("gaus","QL");
	hist->Draw(opt);
	cout << hist->Integral() << "  " << hist->GetRMS() << endl;
	ccdReso.push_back( hist->GetRMS() );
	//ccdReso.push_back( hist->GetFunction("gaus")->GetParameter(2) );
	delete bias; bias=0;
	delete dd; 
}

TGraph *gCcdReso;
void plotSpread() {
	ccdResolution("data/canary_0057.root"); expot.push_back(0);
	ccdResolution("data/canary_0056.root","same"); expot.push_back(0.5);
	ccdResolution("data/canary_0055.root","same"); expot.push_back(1);
	ccdResolution("data/canary_0054.root","same"); expot.push_back(2);
	ccdResolution("data/canary_0053.root","same"); expot.push_back(5);
	ccdResolution("data/canary_0052.root","same"); expot.push_back(10);
	ccdResolution("data/canary_0051.root","same"); expot.push_back(20);
	gCcdReso = new TGraph( expot.size(), &expot[0], &ccdReso[0]);
	
	gCcdReso->Draw("AP");
}





//////////////////////////////////////////////////////////////////////////
//
//	LENGTH CALIBRATION
//
//////////////////////////////////////////////////////////////////////////
float lengthCalibration=0.8547; //  mm/bin 
void length() {
  float y[]={28, 52, 75, 98, 122};
  float x[]={0,   2,  4,  6,   8};
  TGraph *g=new TGraph(sizeof(x)/sizeof(float), x, y);
  g->Draw("AP");
  g->Fit("pol1");
  lengthCalibration=1./g->GetFunction("pol1")->GetParameter(1)*10; // mm/bin
  cout << "New length calibration constant = " << lengthCalibration << " mm/bin" << endl; 
}



// find vertex position and energy lost in travel from vertex to vewifield
float sourceX=130.8, sourceXErr=0.8;
float sourceY=63.13, sourceYErr=0.06;
void activity(int minrun, int maxrun);
void vertex(int irun) {
  irun++;
  //activity(irun, irun);
  gStyle->SetOptFit(1);
  
  TFile *fin = new TFile("canaryCounting.root");
  TTree *ft=(TTree*)fin->Get("ft");
  TProfile *hprof= new TProfile("hprof","",100, -1, 1);
  //TProfile *hprof= new TProfile("hprof","",100, -0.2, 0.2); // WIMP
  //ft->Draw("b:slope>>hprof","b>10");
  ft->Draw("b:slope>>hprof","fyield>35000&&fyield<50000&&b<2000&&slope<0.2");
  hprof->Fit("pol1");
  sourceX=-hprof->GetFunction("pol1")->GetParameter(1);
  sourceY=hprof->GetFunction("pol1")->GetParameter(0);
  sourceXErr=-hprof->GetFunction("pol1")->GetParError(1);
  sourceYErr=hprof->GetFunction("pol1")->GetParError(0);
  hprof->GetFunction("pol1")->SetParNames("source Y", "source X");
  hprof->SetXTitle("a");
  hprof->SetYTitle("b");
  cout << "Source x ..... " << sourceX << " +/- " << sourceXErr << endl;
  cout << "Source y ..... " << sourceY << " +/- " << sourceYErr << endl;
}





//////////////////////////////////////////////////////////////////////////
//
//	GAIN STUDY
//
//////////////////////////////////////////////////////////////////////////
time_t t0=0;
MaxCamSRIM *he=0;
int xmin=100, xmax=110;
void measureGain(int irun) {

  TString fname= irun<10 ? "data/canary_000" : "data/canary_00";
  fname += irun;
  fname += ".root";
  	
  d= new DmtpcDataset;
  d->openRootFile( fname );
  
  TF1 *fun=new TF1("fun","gaus(0)+[3]");
  fun->SetParLimits(0, 1e-10, 1e5);
  fun->SetParLimits(2, 1e-10, 10);
  
  
  float time=0;
  int nev=d->tree()->GetEntries();
  for (int i=0; i<nev; i++) {
    d->getEvent(i);
    
    if (!t0) t0 =  d->event()->timeStamp()->Convert();
    time=(d->event()->timeStamp()->Convert()-t0)/3600.;
    //d->event()->timeStamp()->Print();
    TH2F *image = (TH2F*)d->event()->ccdData(iCCD)->Clone("image");
    if (bias) image->Add(bias, -1);
    TH1D *hy=image->ProjectionY("hy",xmin,xmax);
    hy->GetXaxis()->SetRange(40,80);
    fun->SetParameters(hy->GetMaximum()-hy->GetMinimum(), hy->GetBinCenter( hy->GetMaximumBin()), 1, hy->GetMinimum());
    hy->Fit("fun", "Q");
    //hy->Draw(); c1->Update(); getchar();
    float y    = fun->GetParameter(0) * sqrt2pi * fun->GetParameter(2) / hy->GetXaxis()->GetBinWidth(1);
    float yErr = y * sqrt( pow(fun->GetParError(0)/fun->GetParameter(0),2) + pow(fun->GetParError(2)/fun->GetParameter(2),2) );
    ((TNtuple*)nt)->Fill(y,	yErr,
			 fun->GetParameter(1),	fun->GetParError(1),
			 fun->GetParameter(2),	fun->GetParError(2),
			 fun->GetParameter(3),	fun->GetParError(3),
			 time, irun, xmin, xmax);
    delete image;
    delete hy;
  }
  delete fun;
  delete d; d=0;
  
}


void gain(TString outfilename="canary.root") {

  TFile cfile(outfilename,"RECREATE");
  nt= new TNtuple("nt","","y:yErr:mean:meanErr:width:widthErr:bkg:bkgErr:time:run:xmin:xmax");
  
  makeBiasFrame("data/canary_0001.root");
  
  for (int i=2; i<=32; i++) {
    measureGain(i);
  }
  
  cfile.cd();
  nt->Write();
  cfile.Write();
}


vector<float> ccdGain, ccdGainErr, timeStamp, timeStampErr;
void calcGain_Segment(int irun) {

  if (!he) { he=new MaxCamSRIM("SRIM_He_in_CF4_100Torr"); he->setPressure(75); }
  float x0 = (sourceX-xmax)*lengthCalibration;
  float de=he->calcEnergyLoss(5350, x0, x0+(xmax-xmin+1)*lengthCalibration);
  cout << "Visible energy loss = " << de << "  x0="<<x0<<"  dx="<< (xmax-xmin+1)*lengthCalibration << endl;

  TFile *f= new TFile("canary_gain_run0002-0032_V1-3.root");
  TString cut="width>0.8&&width<1.2&&abs(mean-62.5)<5&&abs(run-";
  cut+=irun;
  cut+=")<0.1";
  TNtuple *nt=(TNtuple*)f->Get("nt");
  TH1F *hy= new TH1F("hy","",100,0, 700e3);
  nt->Draw("y>>hy",cut);
  hy->GetXaxis()->SetRange( hy->GetMaximumBin()-4, hy->GetMaximumBin()+4);
  hy->Fit("gaus","Q");
  ccdGain.push_back(hy->GetFunction("gaus")->GetParameter(1) / de);
  ccdGainErr.push_back(hy->GetFunction("gaus")->GetParError(1) / de);
  
  TH1F *htime=new TH1F("htime","",720,0,360); 
  nt->Draw("time>>htime",cut);
  timeStamp.push_back(htime->GetMean());
  timeStampErr.push_back(htime->GetRMS());
  cout << htime->GetMean()<<"  "<< hy->GetFunction("gaus")->GetParameter(1) << endl;
}

TGraph* plotGain_Segment() {
	for (int irun=2; irun<=32; irun++) calcGain_Segment(irun);
      
	TGraphErrors *g=new TGraphErrors(ccdGain.size(), &timeStamp[0], &ccdGain[0], &timeStampErr[0], &ccdGainErr[0] );
	//g->Fit("pol2");
	g->Draw("AP");
	g->GetHistogram()->SetXTitle("Time (h) ");
	g->GetHistogram()->SetYTitle("Gain (ADU/keV)");
	return g;
}

TString ecalibString="(197.5)";
//TString ecalibString="(197.5-0.887*time+0.00394*time*time)";
//TString ecalibString="(197.5-0.887*(time-1192.758e3)+0.00394*(time-1192.758e3)**2)";

TProfile* plotGain_Full() {

  float xmax=125;
  float xmin=1;

  if (!he) { he=new MaxCamSRIM("SRIM_He_in_CF4_100Torr"); he->setPressure(75); }
  float x0 = (sourceX-xmax)*lengthCalibration;
  float de=he->calcEnergyLoss(5350, x0, x0+(xmax-xmin+1)*lengthCalibration);
  cout << "EnergyLoss = " << de << "  x0="<<x0<<"  dx="<< (xmax-xmin+1)*lengthCalibration << endl;
  
  TFile *f = new TFile("canaryCounting_run0001-0032_V1-2.root");
  TTree *ft=(TTree*)f->Get("ft");
  TProfile *gft=new TProfile("gft","",60,0, 120);
  TString comm="fyield/"; comm+=de; comm+=":time-1192.758e3>>gft";
  ft->Draw( comm, "width<1.3&&width>0.8&&span>70&&yield<7e5&&fyield>1e5","prof");
  //TGraph *gft = (TGraph*)gPad->GetPrimitive("Graph")->Clone("gft");
  //gft->Fit("pol2");

  return gft;
}


void plotGain_All() {
  TProfile *hp=plotGain_Full();
  hp->SetMarkerStyle(25);
  hp->SetXTitle("Time (h)");
  hp->SetYTitle("Gain (ADU/keV)");
  TGraph *gp=plotGain_Segment();
  gp->SetMaximum(250);
  gp->SetMinimum(0);
  hp->Draw("same");
  
  TLegend *leg= new TLegend(0.5, 0.4, 0.8, 0.5);
  leg->AddEntry(gp, "Track segment", "P");
  leg->AddEntry(hp, "Full track", "P");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw();

}

void plotMulti_Full() {

  float xmax=125;
  float xmin=1;

  if (!he) { he=new MaxCamSRIM("SRIM_He_in_CF4_100Torr"); he->setPressure(75); }
  float x0 = (sourceX-xmax)*lengthCalibration;
  float de=he->calcEnergyLoss(5350, x0, x0+(xmax-xmin+1)*lengthCalibration);
  cout << "EnergyLoss = " << de << "  x0="<<x0<<"  dx="<< (xmax-xmin+1)*lengthCalibration << endl;
  
  //TFile *f = new TFile("canaryCounting_run0001-0032_V1-2.root");
  TFile *f = new TFile("canaryCounting_run0058-0077_V1-3.root");
  TTree *ft=(TTree*)f->Get("ft");
  TH2F *gft=new TH2F("gft","",36,0.5, 36.5, 100, 0, 500);
  TString comm="fyield/"; comm+=de; comm+=":run>>gft";
  TString cut="width<1.3&&width>0.8&&span>70&&yield<7e5&&fyield>1e5";
  ft->Draw( comm, cut , "box");

  TH1D *hmulti = gft->ProjectionX();
  hmulti->SetXTitle("Run number");
  hmulti->SetYTitle("Tracks");
  hmulti->Fit("pol1");
  hmulti->Draw("E");
}





TGraph* plotEnergyRange() {
  
  TFile *f = new TFile("canaryCounting_run0001-0032_V1-4.root");
  TTree *ft=(TTree*)f->Get("ft");
  TString comm="span*"; comm+=lengthCalibration; comm+=":fyield/"; comm+=ecalibString; comm +=">>her";
  TH2F *her=new TH2F("her","",100, 0, 4000, 100, 0, 120);
  ft->Draw( comm, "width<1.3&&width>0.8&&span>0&&yield<7e5&&fyield>-1e5&&abs(slope)<0.3");
  her->Draw("box");
  her->SetXTitle("Energy (keV)");
  her->SetYTitle("Range (mm)");

  if (!he) { he=new MaxCamSRIM("SRIM_He_in_CF4_100Torr"); he->setPressure(75); }
  he->getRangeVsEnergy()->Draw("L");

  return 0;
}





//////////////////////////////////////////////////////////////////////////
//
//	TRACK COUNTING
//
//////////////////////////////////////////////////////////////////////////

struct EVENT_DATA {
  float time;
  int   run;
  int   nevents;
  int   sparkEvents;
  int   alphaEvents;
  float eventExposure;
  float volume;
  float pressure;
} tdata;

struct TRACK_DATA {
  int run;
  int event;
  int edge_min_x; // pixels
  int edge_min_y; // pixels
  int edge_max_x; // pixels
  int edge_max_y; // pixels
  int ntrack; // total tracks (>=0) in this event; -1=spark
  float span; // in pixels
  float yield; // in ADU
  float slope; // fitted slope
  float time; // time in h
  float fyield; // fitted yield
  float width; // width of residuals
  float intercept; // fitted intercept y(x=0)
  float eMom[4]; // energy moments for bins (mean, rms, skew, curt)
  float rMom[4]; // range moments for bins
  float maxPix; // maximum energy in pixel  
  void print() {
    cout << run << ","<< event << ": (" 
	 <<  edge_min_x<<","<<edge_min_y<<"),("
	 <<  edge_max_x<<","<<edge_max_y<<")  R=" << span
	 << "  y="<<yield 
	 << "  yfit=" << fyield
	 << "  w=" << width
	 << "  a="<<slope
	 << "  b="<<intercept
	 << "  num.tracks="<< ntrack
	 << "  T="<<time <<" h"  
	 << endl;
  }
} fitData;

TH1F *hyield=0;

void findTracks(int irun, float threshold=900, bool doPlot=false) {

  //threshold=30;//WIMP

  float ceiling=14000;  // 16384

  // Detector setup
  tdata.time=-1;
  tdata.run=irun;
  tdata.alphaEvents=0;
  tdata.sparkEvents=0;
  tdata.eventExposure=5; // sec
  tdata.volume=1; // lit
  tdata.pressure=75; // Torr
  
  TDatime time0(2009, 5, 19, 16, 0, 0);

  // Input name
  TString fname="data/canary_";
  //TString fname="data/dmtpc_run0"; // WIMP
  if (irun<1000) fname+="0";
  if (irun<100) fname+="0";
  if (irun<10) fname+="0";
  fname += irun;
  fname += ".root";
  d= new DmtpcDataset;
  d->openRootFile( fname );
  TH2F *bias=d->getBiasFrame(iCCD+1);
  
  // Viewing pad
  TCanvas *c=new TCanvas("c","",900, 600);
  c->Divide(3,2);

  // Spark counting
  float sparkThreshold=2000;
  //float sparkThreshold=1065; // WIMP

  // alpha yield
  TF1 *fun=new TF1("fun","gaus(0)+[3]");

  // Loop over events
  tdata.nevents=d->tree()->GetEntries();
  for (int i=0; i<tdata.nevents; i++) {
    d->getEvent(i);
    TH2F *image = (TH2F*)d->event()->ccdData(iCCD)->Clone("image");

    // get time stamp
    if (tdata.time<0) tdata.time = (d->event()->timeStamp()->Convert()-time0.Convert())/3600.;

    // is spark?
    bool isSpark = d->event()->ccdData(iCCD)->Integral()/d->event()->ccdData(iCCD)->GetEntries()>sparkThreshold ? true : false; 
    if (isSpark) { cout << "*"; tdata.sparkEvents++; }
    else cout << ".";
    cout << flush;

    // raw image
    c->cd(1);
    if (bias) image->Add(bias, -1);
    if (doPlot) image->DrawCopy("colz");
    TH2F *rawimage=(TH2F*)image->Clone("rawimage");   
 
    // clean up image
    MaxCamImageTools::applyThreshold(image, threshold);
    MaxCamImageTools::applyCeiling(image, ceiling);
    MaxCamImageTools::killLonePixels(image, threshold, 4);
    c->cd(2);
    if (doPlot) image->DrawCopy("colz"); 

  
    // track fit
    if (!isSpark) {
      c->cd(3);
      MaxCamTrack *trfit = new MaxCamTrack( image, false);
      trfit->setThreshold( threshold );
      trfit->setCleaningSpan(5);
      trfit->makeTracks();
      if (doPlot) {
	image->DrawCopy("colz");
	for (int it=0; it<trfit->nTracks(); it++) {
	  trfit->getTrack(it)->DrawCopy("same");
	}
      }
      if ( trfit->nTracks()>0 ) {
	tdata.alphaEvents++;
	
	// cannot handle multiple-alpha tracks so
	// take all info from 1st track or whole image
	int imax, jmax;
	MaxCamImageTools::countSegmentLength2D(image, threshold, imax, jmax);
	
	fitData.run=irun;
	fitData.event=i;
	if (MaxCamImageTools::distanceToImageEdge(image, imax) > 
	    MaxCamImageTools::distanceToImageEdge(image, jmax)) {
	  int tmp=imax;
	  imax=jmax;
	  jmax=tmp;
	}
	fitData.edge_min_x=imax%(image->GetNbinsX()+2);
	fitData.edge_min_y=imax/(image->GetNbinsX()+2);
	fitData.edge_max_x=jmax%(image->GetNbinsX()+2);
	fitData.edge_max_y=jmax/(image->GetNbinsX()+2);
	fitData.ntrack=trfit->nTracks();
	fitData.span=MaxCamImageTools::calcPixelDistance(image, imax, jmax);
	fitData.yield=image->Integral();
	fitData.slope=trfit->getTrack(0)->GetParameter(1);
	fitData.intercept=trfit->getTrack(0)->Eval(0);
	fitData.time=tdata.time+i*tdata.eventExposure/3600.;
	fitData.maxPix=image->GetMaximum();

	c->cd(4);
	if (1||!hyield) hyield=MaxCamImageTools::createYieldHisto( image, threshold, ceiling );  
	else hyield->Add( MaxCamImageTools::createYieldHisto( image, threshold, ceiling ) );  
       //TH1F *hyield=MaxCamImageTools::createYieldHisto( image, threshold, ceiling );
	if (doPlot) hyield->DrawCopy();
	fitData.eMom[0]=hyield->GetMean();
	fitData.eMom[1]=hyield->GetRMS();
	fitData.eMom[2]=hyield->GetSkewness();
	fitData.eMom[3]=hyield->GetKurtosis();

	// residuals
	c->cd(5);
	trfit->setImage(rawimage);
	trfit->setThreshold( threshold*0.25 );
	TH1F *hres=trfit->makeResiduals(0, "th", 10);
	//TH1F *hres=trfit->makeResiduals(0, "th", 30); // WIMP
	fun->SetParameters(hres->GetMaximum(), 0, 1, hres->GetMinimum());
	hres->Fit("fun", doPlot ? "" : "Q0" );
	TF1 *ffun=hres->GetFunction("fun");
	fitData.width=fabs(ffun->GetParameter(2));
	fitData.fyield=ffun->GetParameter(0)*fitData.width*sqrt2pi/hres->GetXaxis()->GetBinWidth(1);

	fitData.print();
	ft->Fill();


	if (doPlot) { c->Update(); getchar(); }

	delete hyield;
	delete hres;
	delete image;
	delete trfit;
      }
      else {
	fitData.ntrack = 0;
	ft->Fill();

	delete trfit;
	delete rawimage;
      }

    } else {
      fitData.ntrack = -1;
      ft->Fill();      
    }


    //c->Update(); getchar();
  }
  delete d; d=0;
  delete fun;

  float activity = tdata.alphaEvents/( tdata.eventExposure * tdata.nevents * tdata.volume );
  float activityErr = sqrt(tdata.alphaEvents)/( tdata.eventExposure * tdata.nevents * tdata.volume );

  cout << "Run ................. " << tdata.run << endl;
  cout << "Total exposure ...... " << tdata.eventExposure*tdata.nevents/3600 << " h" << endl;
  cout << "Total volume ........ " << tdata.volume << " lit"<< endl;
  cout << "Total mass .......... " << MaxCamSRIM::density(75)*tdata.volume*1e3 << " g" <<endl;
  cout << "Alpha events ........ " << tdata.alphaEvents << endl;
  cout << "Sparks .............. " << tdata.sparkEvents << endl;
  cout << "Total activity ...... " << activity << " +- " << activityErr << " Bq/lit"<< endl;
  cout << "               ...... " << activity/37e9*1e12 << " +- " << activityErr/37e9*1e12 << " pCi/lit"<<endl;

  nt->Fill();
}



void activity(int minrun, int maxrun) {

  // output file
  TFile outfile("canaryCounting.root","UPDATE");
  //TFile outfile("test.root","RECREATE");
  nt = (TTree*)outfile.Get("nt");
  ft = (TTree*)outfile.Get("ft");
  if (!nt) {
	nt = new TTree("nt","");
  	nt->Branch("data", &tdata, "time/F:run/I:nevents:sparkEvents:alphaEvents:exposure/F:volume:pressure");
  }
  else {
    nt->SetBranchAddress("data",&tdata);
  }
  if (!ft) {
	ft = new TTree("ft","");
  	ft->Branch("fit", &fitData, "run/I:event:minx:miny:maxx:maxy:ntrack:span/F:yield:slope:time:fyield:width:b:eMom[4]:rMom[4]:maxPix");
  }
  else ft->SetBranchAddress("fit",&fitData);

  // loop over runs
  for (int irun=minrun; irun<=maxrun; irun++) findTracks(irun);

  // save
  outfile.cd();
  nt->Write();
  ft->Write();
  outfile.Write();
}



void plotActivity_Rough(TString what="alpha") {
  //TFile infile("canaryCounting_run0058-0077_V1-1.root");
  TFile infile("canaryCounting_tmp.root");
    TTree *nt=(TTree*)infile.Get("nt");
    nt->SetBranchAddress("data",&tdata);
    vector<float> time, timeErr, spark, sparkErr, alpha, alphaErr;
    for (int i=0; i<nt->GetEntries(); i++) {
      nt->GetEvent(i);
      time.push_back( tdata.time );
      timeErr.push_back( 0 );
      spark.push_back( tdata.sparkEvents );
      sparkErr.push_back( sqrt(tdata.sparkEvents) );
      alpha.push_back( tdata.alphaEvents );
      alphaErr.push_back( sqrt(tdata.alphaEvents) );
      cout << tdata.run 
	   << ": time=" << tdata.time 
	   << "h  alpha=" << tdata.alphaEvents 
	   << "   spark="<<tdata.sparkEvents  << endl;
    }
    TGraphErrors *g=0;
    bool doFit=true;
    if (what=="alpha") g = new TGraphErrors( time.size(), &time[0], &alpha[0], &timeErr[0], &alphaErr[0] );
    else if (what=="spark")  g = new TGraphErrors( time.size(), &time[0], &spark[0], &timeErr[0], &sparkErr[0] );
    else if (what.Contains("spark") && what.Contains("alpha")) {
      g = new TGraphErrors( alpha.size(), &alpha[0], &spark[0], &alphaErr[0], &sparkErr[0] );
      doFit=false;
    }
    g->Draw("AP");
    if (doFit) g->Fit("pol1");
}


double nonzeroPoisson(double n, double mu) {  return n-mu/(1-exp(-mu)); }
double nonzeroPoissonDeriv(double mu) {  return -1./(1-exp(-mu))*(1-mu*exp(-mu)/(1-exp(-mu))  ); }
double nonzeroPoissonSolver(double n) {
  
  double oldmu=-1e10, mu=n;
  while (1) {
    mu -= nonzeroPoisson(n,mu)/nonzeroPoissonDeriv(mu);
    //cout << oldmu << "  " << mu << endl;
    if ( fabs(oldmu-mu)<1e-4 ) break;
    oldmu=mu;
  }

  return mu;
}




int measureLifetime(TTree *ft, TString cut, int run, float &mu, float &muErr, float &time, float &timeErr) {

  TF1 *fpoisson=new TF1("fpoisson","[0]*TMath::Poisson(x,[1])",0,10);
  TH1F *hnt=new TH1F("hnt","", 10, 0, 10);
  cut += "&&run==";
  cut += run;

  ft->Draw("ntrack>>hnt",cut);
  int nentr=hnt->GetEntries();
  if (!nentr) return 0;
  fpoisson->SetParameters( hnt->GetMaximum(), hnt->GetMean());
  //hnt->Fit("fpoisson","Q");
  //c1->Update();getchar(); 
  mu=hnt->GetMean();//GetFunction("fpoisson")->GetParameter(1);
  muErr=hnt->GetRMS()/sqrt(hnt->GetEntries());//GetFunction("fpoisson")->GetParError(1);

  TH1F *ht=new TH1F("ht","",6000,0, 2000);
  ft->Draw("time>>ht",cut);
  ht->Fit("gaus","Q");
  time= ht->GetMean();//GetFunction("gaus")->GetParameter(1);
  timeErr = ht->GetRMS();//GetFunction("gaus")->GetParError(1);

  delete fpoisson; delete hnt; delete ht;
  return nentr;
}

double fHL(double *x, double *par) {
  
  double t=x[0];
  double fun=0;
  fun += par[0];

  if (t>505 && t<688) {
    fun += par[1]*exp( par[2]*(t-500) );
    fun += par[3]*exp( par[4]*(t-500) );
  }
  
  if (t>890) fun += par[5];
	
	if (t>1030 && t<1127) {
		fun += par[6]*exp( par[7]*(t-1000) );
		fun += par[8]*exp( par[9]*(t-1000) );
	}
	
	if (t>1198) fun += par[10];

  return fun;
}




void plotActivity_Hough(TString what) {


  TChain *ft=new TChain("ft");
  ft->Add("canaryCounting_59-70_Baseline.root");
  ft->Add("canaryCounting_86-90_SS_Background.root");
  ft->Add("canaryCounting_93-114_SS_Radon2.root");
  ft->Add("canaryCounting_115-122_SS_Background2.root");
  ft->Add("canaryCounting_126-134_Cu_Background.root");
  ft->Add("canaryCounting_135-143_Cu_Radon2.root");
  ft->Add("canaryCounting_144-150_Cu_NoRadon.root");
	
  TString cut="width<1.3&&width>0.6&&span>10&&fyield>10&&fyield<1e6";


  if (what.Contains("energy") && what.Contains("range")) {
    TString comm="span*"; comm+=lengthCalibration; comm+=":fyield/"; comm+=ecalibString;
    cut += "&&ntrack==1";
    ft->Draw( comm, cut);
    TGraph *g = (TGraph*)gPad->GetPrimitive("Graph")->Clone("g");
    gPad->Clear();
    g->Draw("AP");
    g->GetHistogram()->SetXTitle("Energy (keV)");
    g->GetHistogram()->SetYTitle("Range (mm)");    
    
    if (!he) { he=new MaxCamSRIM("SRIM_He_in_CF4_100Torr"); he->setPressure(75); }
    he->getRangeVsEnergy()->Draw("L");
  }

  else if ( what.Contains("multi") && what.Contains("run") ) {

    TH2F *gft=new TH2F("gft","",  21, 57.5, 78.5, 100, 0, 5e5);
    TString comm="fyield:run>>gft";
    cut += "&&ntrack==1";
    ft->Draw( comm, cut , "goff");
    TH1D *hmulti = gft->ProjectionX();
    hmulti->SetXTitle("Run number");
    hmulti->SetYTitle("Tracks (1/5h)");
    //hmulti->Fit("pol1");
    hmulti->SetMarkerStyle(24);
    hmulti->DrawCopy("E");

    cut += "&&minx>10&&minx<115&&miny>10&&miny<115";
    ft->Draw( comm, cut , "goff");
    hmulti = gft->ProjectionX();
    hmulti->Draw("same E");
  }
 
  else if (what=="yield") {
    TH2F *gft=new TH2F("gft","",  10, 0, 170, 100, 0, 5e5);
    TString comm="fyield:time>>gft";
    cut += "&&ntrack==1";
    ft->Draw( comm, cut , "box");
    TH1D *hmulti = gft->ProjectionX();
    hmulti->SetXTitle("Run number");
    hmulti->SetYTitle("Tracks");
    hmulti->Fit("pol1");
    hmulti->Draw("E");    
  }

  else if ( what.Contains("multi") && what.Contains("time") ) {

    TProfile *hmulti=new TProfile("hmulti","", 60 , 500, 700);
    TString comm="ntrack:time>>hmulti";
    cut += "&&ntrack>0";
    ft->Draw( comm, cut, "goff");
    hmulti->SetXTitle("Time (h)");
    hmulti->SetYTitle("<Num. track>");
    hmulti->SetMarkerStyle(24);
    hmulti->DrawCopy("E");

    if (what.Contains("corr")) { // correct for 0-track events
      TH1F *hmu=new TH1F("hmu","", hmulti->GetNbinsX(), hmulti->GetXaxis()->GetXmin(), hmulti->GetXaxis()->GetXmax());
      for (int i=1; i<=hmulti->GetNbinsX(); i++) {
	float mu_raw = hmulti->GetBinContent(i);
	if (mu_raw>0) {
	  float mu_corr= nonzeroPoissonSolver(mu_raw);
	  float mu_corr_err=sqrt(mu_corr/hmulti->GetBinEntries(i));
	  //cout << mu_raw << "  " << mu_corr << endl;
	  hmu->SetBinContent(i, mu_corr);
	  hmu->SetBinError(i, mu_corr_err );
	}
      }
      hmu->SetMarkerStyle(20);
      TF1 *fdec=new TF1("fdec","expo(0)", 525, 700);
      fdec->SetParameters(60, -0.1, 2.3, -6e-3);
      hmu->Fit("fdec", "R");
      float HLFactor= 1./24*log(2);
      int HLPar=1;
      cout << "Lifetime ........... " << -HLFactor/fdec->GetParameter(HLPar)/log(2) << " days" << endl;
      cout << "Half-life .......... " << -HLFactor/fdec->GetParameter(HLPar) << " +/- " 
      	   << HLFactor/pow(fdec->GetParameter(1),2)*fdec->GetParError(HLPar)  << " days" << endl;
      hmulti->DrawCopy("E");
      hmu->DrawCopy("same");
    }

  }

  else if (what=="activity") {

    cut += "&&ntrack>=0";
    vector<float> mu, muErr, time, timeErr;
    int run0=59, run1=150;
    for (int irun=run0; irun<=run1; irun++)  {
      if (irun==108) continue;
	  if (irun==144) continue;
      //if (irun>=93 && irun<=96) continue;
      float m,me,t,te;
      int n=measureLifetime(ft, cut, irun, m, me, t, te);
      if (!n) continue;
      mu.push_back(m); muErr.push_back(me); time.push_back(t); timeErr.push_back(te);
      cout << "run="<<irun<<"  act="<< m << "+/-" << me << "   " << t << "+/-" << te << endl;
    }

    TGraphErrors *g=new TGraphErrors(mu.size(), &time[0], &mu[0], &timeErr[0], &muErr[0]);
    g->Draw("AP");
    TF1 *rnfun= new TF1("rnfun", fHL, 0, 2000, 10);
    rnfun->SetParameters(8e-2, 0.24, -7e-3,  0.8, -7e-2, 4e-2 );
    
    g->Fit("rnfun");
    float HLFactor= 1./24*log(2);
    TF1 *fdec=g->GetFunction("rnfun");
    int hlpar=2;
    cout << "Lifetime ........... " << -HLFactor/fdec->GetParameter(hlpar)/log(2) << " days" << endl;
    cout << "Half-life .......... " << -HLFactor/fdec->GetParameter(hlpar) << " +/- " 
	 << HLFactor/pow(fdec->GetParameter(hlpar),2)*fdec->GetParError(hlpar)  << " days" << endl;

  }
 
  

}




void plotVisualCounting() {
  float t[]={5*1000./3600, 5*1500./3600, 5*2000./3600, 5*2500./3600, 5*3000./3600,  5*3500./3600,  
	     13, 13+5*1000./3600, 13+5*2000./3600, 13+5*3500./3600};
  float terr[]={0,0,0,0,0,0,0,0,0};
  float y[]={88./21, 70./25, 54./26, 38./26, 27./22,  38./30,  31./29, 26./26, 24./26, 21./24};
  float yerr[]={21, 25, 26, 26, 22, 30, 29, 26, 26, 24};
  int n=sizeof(y)/sizeof(float);
  for (int i=0; i<n; i++)   yerr[i]=sqrt(y[i]/yerr[i]);
  TGraphErrors *g= new TGraphErrors(n,t,y, terr,yerr);
  g->Draw("AP");
}




//////////////////////////////////////////////////////////////////////////
//
//	WORM REJECTION
//
//////////////////////////////////////////////////////////////////////////



void mlpWorm(Int_t ntrain=100) {
   
   if (!gROOT->GetClass("TMultiLayerPerceptron")) {
      gSystem->Load("libMLP");
   }

   // Prepare inputs
   // The 2 trees are merged into one, and a "type" branch, 
   // equal to 1 for the signal and 0 for the background is added.
   Int_t type;

   TFile *fAlpha = TFile::Open( "canaryCounting_run0002.root" );
   TTree *alphas = (TTree *) fAlpha->Get("ft");
   alphas->SetBranchAddress("fit",  &fitData);

   TFile *fWorm = TFile::Open( "canaryCounting_run0086-0089.root" );
   TTree *worms = (TTree *) fWorm->Get("ft");
   worms->SetBranchAddress("fit",  &fitData);

   TFile *fNN = TFile::Open( "canaryCounting_nntrain.root", "RECREATE" );
   TTree *simu = new TTree("MonteCarlo", "Filtered Monte Carlo Events");
   simu->Branch("fit", &fitData, "run/I:event:minx:miny:maxx:maxy:ntrack:span/F:yield:slope:time:fyield:width:b:eMom[4]:rMom[4]:maxPix");
   simu->Branch("type",   &type,   "type/I");
   fNN->Print();
   
   type = 1; for (int i = 0; i < alphas->GetEntries(); i++) { alphas->GetEntry(i); if (fitData.eMom[0]>0) simu->Fill(); }
   type = 0; for (int i = 0; i < worms->GetEntries(); i++) { worms->GetEntry(i); if (fitData.eMom[0]>0) simu->Fill(); }

/*
   TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron("@eMom[0],@eMom[1],@eMom[2],@eMom[3],@width,@maxPix:12:6:type",
                                                          simu,
							  "Entry$%2",
							  "(Entry$+1)%2");
   mlp->Train(ntrain, "text,graph,update=10");
   mlp->Export("test");


   // Use TMLPAnalyzer to see what it looks for
   TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network analysis");
   mlpa_canvas->Divide(2,2);
   TMLPAnalyzer ana(mlp);
   // Initialisation
   ana.GatherInformations();
   // output to the console
   ana.CheckNetwork();
   mlpa_canvas->cd(1);
   // shows how each variable influences the network
   ana.DrawDInputs();
   mlpa_canvas->cd(2);
   // shows the network structure
   mlp->Draw();
   mlpa_canvas->cd(3);
   // draws the resulting network
   ana.DrawNetwork(0, "type==1", "type==0");

   gStyle->SetMarkerStyle(1);
   TCanvas* varc=new TCanvas("varc","", 900,600); varc->Divide(3,1);
   varc->cd(1); TH2F *fr01=new TH2F("fr01","",100,0,16384, 100,0,8000); fr01->SetLineColor(2);
   alphas->Draw("eMom[1]:eMom[0]>>fr01");  fr01->DrawCopy("box");  
   worms->Draw("eMom[1]:eMom[0]>>fr01","eMom[0]>0","goff"); fr01->SetLineColor(4);  fr01->DrawCopy("same box");  
   
   varc->cd(2); TH2F *fr02=new TH2F("fr02","",100,-5,5, 100,-3,3); fr02->SetLineColor(2);
   alphas->Draw("eMom[2]:eMom[3]>>fr02");  fr02->DrawCopy("box");  
   worms->Draw("eMom[2]:eMom[3]>>fr02","eMom[0]>0","goff"); fr02->SetLineColor(4);  fr02->DrawCopy("same box");  

   varc->cd(3); TH2F *fr03=new TH2F("fr03","",100,0,16384, 100,0,10); fr03->SetLineColor(2);
   alphas->Draw("width:maxPix>>fr03");  fr03->DrawCopy("box");  
   worms->Draw("width:maxPix>>fr03","eMom[0]>0","goff"); fr03->SetLineColor(4);  fr03->DrawCopy("same box");     
*/
}



//////////////////////////////////////////////////////////////////////////
//
//	RADON DECAY CHAINS
//
//////////////////////////////////////////////////////////////////////////


void Rn_Chains() {
  float day=1;
  float hour=1./24;
  float min=hour/60.;
  float sec=min/60.;

  TH2F *frame=new TH2F("frame","", 200, 3*sec, 10*day,  200, 0, 20);
  frame->Draw();

  // Rn222 chain
  float alphaEnergy[] = { 5.5,           6,        0,       0,    7.7};
  float alphaHL[]     = { 3.8*day, 3.1*min, 26.8*min,  20*min,      0.2e-3*sec};
  vector<float> Rn222_Energy;
  vector<float> Rn222_Time;
  int n=sizeof(alphaHL)/sizeof(float);
  float totalTime=0;
  for (int i=0; i<n; i++) {
    totalTime+=alphaHL[i];
    if ( alphaEnergy[i]>0 ) {
      Rn222_Time.push_back( totalTime );
      Rn222_Energy.push_back( alphaEnergy[i] );
    }
  }
  TGraph *g222=new TGraph( Rn222_Time.size(), &Rn222_Time[0], &Rn222_Energy[0]);
  g222->Draw("PL");


  // Rn220 chain
  float alphaEnergy220[] = { 10.6,           1.78,          0,          0,           13.3};
  float alphaHL220[]     = { 55.6*sec,  0.145*sec, 10.64*hour,  60.55*min,         45*sec};
  vector<float> Rn220_Energy;
  vector<float> Rn220_Time;
  n=sizeof(alphaHL220)/sizeof(float);
  totalTime=0;
  for (int i=0; i<n; i++) {
    totalTime+=alphaHL220[i];
    if ( alphaEnergy220[i]>0 ) {
      Rn220_Time.push_back( totalTime );
      Rn220_Energy.push_back( alphaEnergy220[i] );
    }
  }
  TGraph *g220=new TGraph( Rn220_Time.size(), &Rn220_Time[0], &Rn220_Energy[0]);
  g220->Draw("PL");



}
