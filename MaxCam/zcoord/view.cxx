//   root[0] .L view.cxx++
//   root[1] init(filename)
//   root[0] next()
//

#include "TString.h"
#include "TCanvas.h"
#include "../DmtpcDataset.hh"
#include "../MaxCamWaveformTools.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "../ScopeDataInfo.hh"

#include <iostream>
using std::cout;
using std::endl;

DmtpcDataset *d;
TH2F *image[2], *bias[2];
int iev, nEntries;
TCanvas *c_ccd=0, *c_scope=0;
bool firstPass;
// # of channels on the scope.  
// this matters because the wfs from all scope channels 
// are returned in a 1d array
Int_t nChan = 2;  

const int maxWF=3;

void incr();
void decr();
void setiev(int n);
void drawScopeTraces();
void next();
void prev();

void init(TString fname) {
  c_scope = new TCanvas("c_scope","",400,800); 
  c_scope->Divide(1,maxWF);

  d = new DmtpcDataset;
  d->openRootFile(fname);

  nEntries = d->tree()->GetEntries();
  cout << "d->tree()->GetEntries() = " << d->tree()->GetEntries() << endl;

  iev=0;
  firstPass = true;
}


void next() {
  if (!firstPass) incr();
  d->getEvent(iev);
  drawScopeTraces();
  firstPass = false;
}

void prev() {
  decr();
  d->getEvent(iev);
  drawScopeTraces();
}

void update() {
  drawScopeTraces();
}

void drawScopeTraces() {
  cout << "iev = " << iev << "/" << (nEntries-1) << endl;
  int nwf=d->event()->scopeData()->GetEntries()/nChan;
  int iwf=0;
  Int_t nWfOmitted = 0;
  for (int ii=0; ii<nwf; ii++) {
    if (ii>(maxWF-1)) { 
      nWfOmitted++;
      continue;
    }
    c_scope->cd(ii+1); 
    d->event()->scopeData(ii)->Draw();
    //c_scope->cd(2*iwf+1); d->event()->scopeData(i)->Draw();
    //c_scope->cd(2*iwf+2); d->event()->scopeData(i+nwf)->Draw();
    MaxCamWaveformTools wfA(d->event()->scopeData(ii)); wfA.print();
    //MaxCamWaveformTools wfB(d->event()->scopeData(i+nwf)); wfB.print();
    iwf++;
  }
  cout << "total wf="<<nwf << " plotted="<<iwf<<endl;
  if (nWfOmitted > 0) {
    cout << "******* " << nWfOmitted << " waveform(s) not drawn ********   " << endl;
  }
}


void decr() {
  if (iev > 0) iev--;
}

void incr() {
  if (iev < (nEntries-1)) iev++;
}
void setiev(int n) {
  iev = n;
}

void loop(int n) {
  for (int i=0; i<n; i++) { 
    cout << "event="<<iev << endl;
    next(); 
  }
}


void nTrigSpectrum(TString fname) {
  // make a spectrum of the number of triggers per "exposure"
  // also make a plot of baseline vs time

  TCanvas *ctrg = new TCanvas("trgspec","",400,400); 

  d = new DmtpcDataset;
  d->openRootFile(fname);

  nEntries = d->tree()->GetEntries();
  cout << "d->tree()->GetEntries() = " << d->tree()->GetEntries() << endl;

  // FIX hist range...
  Int_t trgMin = 0;
  Int_t trgMax = 20;
  TH1F *ntrgh = new TH1F("ntrgh", "", trgMax-trgMin+1, trgMin, trgMax); 

  // loop over all events
  Int_t nTrigInEvent;
  Int_t nTrigTotal = 0; 
  for (Int_t iev=0; iev<nEntries; iev++) {
    //cout << "event: " << iev << endl;
    d->getEvent(iev);
    // get number of triggers
    nTrigInEvent = d->event()->scopeData()->GetEntries()/nChan;
    for (Int_t iwf=0; iwf<nTrigInEvent; iwf++) {
      //cout << "   iwf = " << iwf << endl;
      // loop over wfs and get pulse heights -- push into a TH1F
      MaxCamWaveformTools wfA(d->event()->scopeData(iwf)); 
      //wfA.print();
      //cout << "        " << wfA.getPeakHeight() << endl;
      //mca1->Fill(1e3*wfA.getPeakHeight());
      ntrgh->Fill(1e3*wfA.getTroughDepth());
      nTrigTotal++;
    }
  }

  ntrgh->Draw();
  ntrgh->GetXaxis()->SetTitle("Peak Voltage (mV)");
  ntrgh->GetYaxis()->SetTitle("#");
  ctrg->Update();
  cout << "there were " << nTrigTotal << " triggers" << endl;

  TFile outfile("ntrgHist.root", "RECREATE");
  ntrgh->Write();
  outfile.Close();
}

void mca(TString fname) {

  c_scope = new TCanvas("scopemca","mca",400,400); 

  d = new DmtpcDataset;
  d->openRootFile(fname);

  nEntries = d->tree()->GetEntries();
  cout << "d->tree()->GetEntries() = " << d->tree()->GetEntries() << endl;

  // FIX hist range...
  //TH1F *mca1 = new TH1F("mca1", "Scope MCA", 100, -20., 20.); // millivolts
  TH1F *mca1 = new TH1F("mca1", "Scope MCA", 1000, -100., 100.); // millivolts

  // loop over all events
  Int_t nTrigInEvent;
  Int_t nWfInEvent;
  Int_t nTrigTotal = 0; 
  for (Int_t iev=0; iev<nEntries; iev++) {
    //cout << "event: " << iev << endl;
    d->getEvent(iev);
    // get number of triggers
    nTrigInEvent = d->event()->scopeData()->GetEntries();
    nWfInEvent   = nTrigInEvent/nChan;
    for (Int_t iwf=0; iwf<nWfInEvent; iwf++) {
      //cout << "   iwf = " << iwf << endl;
      // loop over wfs and get pulse heights -- push into a TH1F
      MaxCamWaveformTools wfA(d->event()->scopeData(iwf)); 
      //wfA.print();
      //cout << "        " << wfA.getPeakHeight() << endl;
      //mca1->Fill(1e3*wfA.getPeakHeight());
      mca1->Fill(1e3*wfA.getTroughDepth());
      nTrigTotal++;
    }
  }

  mca1->Draw();
  mca1->GetXaxis()->SetTitle("Peak Voltage (mV)");
  mca1->GetYaxis()->SetTitle("#");
  cout << "there were " << nTrigTotal << " triggers" << endl;

  TFile outfile("mcaHist.root", "RECREATE");
  mca1->Write();
  outfile.Close();
}


void allmcas() {
  TString fileroot("sequence_2009-11-12T18:04:29_");

  TCanvas *c1 = new TCanvas();
  TString firstfile(fileroot);
  firstfile += "0_mca.root";
  TFile *ff = new TFile(firstfile);
  TH1F* htemp = (TH1F*)ff->Get("mca");
  htemp->SetLineStyle(1);
  htemp->Draw();

  for (Int_t ii=1; ii<10; ii++) {
    TString filename(fileroot);
    filename += ii;
    filename += "_mca.root";
    cout << filename << endl;
    TFile *ff = new TFile(filename);
    TH1F* htemp = (TH1F*)ff->Get("mca");
    htemp->SetLineStyle(ii+1);
    htemp->Draw("SAME");
    //((TH1F*)ff->Get("mca"))->Draw("SAME");
  }

  c1->Update();
  c1->Print("junk.pdf");
}
