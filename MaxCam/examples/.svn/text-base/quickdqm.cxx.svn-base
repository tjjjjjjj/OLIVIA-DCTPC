#include "../DmtpcDataset.hh"
#include "../MaxCamChannel.hh"
#include "../DmtpcEvent.hh"
#include "../MaxCamImageTools.hh"

#include "TH2F.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TDatime.h"
#include "TFile.h"
#include "TCanvas.h"

#include "time.h"
#include "iostream"
using std::cout;
using std::endl;


DmtpcDataset *d;
TH2F *image, *bias;
TF1 *fgauss;
int iev;
int iccd=0;
TNtuple *nt;
time_t T0, Tprev;
TFile *of;

void getBiasFrame() {
  DmtpcDataset bd;
  bd.openRootFile("data/dmtpc_run00321.root");
  bias=(TH2F*)bd.getBiasFrame(1)->Clone("bias");
}


void init(int irun) {

  //getBiasFrame();

  d= new DmtpcDataset;
  TString fname="dmtpc_run";
  if (irun<100) fname+="000";
  else if (irun<1000) fname+="00";
  else if (irun<10000) fname+="0";
  fname += irun;
  TString fnameRead="data/";
  fnameRead += fname;
  fnameRead += ".root";

  d->openRootFile(fnameRead);
  iev=0;
  d->getEvent(0);
  T0=d->event()->timeStamp()->Convert();
  Tprev=T0;
  fgauss=new TF1("fgauss","gaus");
  bias=d->getBiasFrame(1);

  TString fnameWrite = fname;
  fnameWrite += "_DQM_Ver001.root";
  of = new TFile(fnameWrite,"CREATE");
  nt=new TNtuple("nt","","index:mean:width:skewness:kurtosis:n3sigma:dt:t:anode:drift:press");
}


void save() { of->Write(); }



void analyze() {
  image=d->event()->ccdData(iccd);
  if (bias) image->Add(bias,-1);
  //image->SetMinimum(-50); image->SetMaximum(150); image->Draw("colz");
  TH1F *yield=MaxCamImageTools::createYieldHisto(image,-50, 150);
  fgauss->SetParameters( yield->GetMaximum(), yield->GetBinCenter(yield->GetMaximumBin()), 6.7);
  yield->Fit("fgauss","Q");
  float mean = fgauss->GetParameter(1);
  float width = fgauss->GetParameter(2);
  float th=mean+3*width;
  float skewness=yield->GetSkewness();
  float kurtosis=yield->GetKurtosis();
  //yield->DrawCopy();
  delete yield;
  yield=MaxCamImageTools::createYieldHisto(image,th, 150); 
  float n3sigma = yield->Integral();
  //cout << "threshold=" << th << "   n="<<n3sigma<<endl;
  delete yield;
  time_t T=d->event()->timeStamp()->Convert();
  float	anode=d->event()->experimentConfig("anodeHV")->currentValue;
  float	drift=d->event()->experimentConfig("driftHV")->currentValue;
  float	press=d->event()->experimentConfig("pressure")->currentValue;
  nt->Fill(iev, mean, width, n3sigma, skewness, kurtosis, T-Tprev, T, anode, drift, press);
  Tprev=T;
}

void next() {
  d->getEvent(iev++);
  analyze();
}

void loop(int n=0) {
  if (!n) n=d->chain()->GetEntries()-iev;
  for (int i=0; i<n; i++) { 
    if (i%100==0) cout << "event="<<iev << endl;
    next(); 
    //c1->Update();
  }
}
