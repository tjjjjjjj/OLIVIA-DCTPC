//   root[0] .L analyze.cxx++
//

#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "../DmtpcDataset.hh"
#include "../MaxCamWaveformTools.hh"
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

const int nfiles = 6;
TString files[nfiles] = {"zcoord_X_000010.root", 
			 "zcoord_X_000005.root",
			 "zcoord_X_000006.root",
			 "zcoord_X_000007.root",
			 "zcoord_X_000008.root",
			 "zcoord_X_000009.root"};
Float_t means[nfiles];
Float_t pressures[nfiles] = { 5.15,  11.6,  16.6,  22.7,  35.4,  50.2}; // torr
Float_t vfitMin[nfiles]   = {-20.0, -20.0, -20.0, -20.0, -20.0, -20.0};
Float_t vfitMax[nfiles]   = { -2.0, -10.0, -10.0,  -9.5,  -6.0,  -2.0};
Int_t mcaNumber = 0;

const int maxWF=3;

void init(TString fname);
void heightVsPressure();
void mca(TString fname);

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

void heightVsPressure() {

  TCanvas *c_pres = new TCanvas("c_pres", "", 400, 800);
  c_pres->Divide(2, (nfiles/2)+1);

  for (int ii=0; ii<nfiles; ii++) {
    cout << files[ii] << endl;
    c_pres->cd(ii+1);
    mcaNumber = ii;
    mca(files[ii]);
  }

  for (Int_t ii=0; ii<nfiles; ii++) {
    cout << pressures[ii] << "   " << means[ii] << endl;
  }

  TGraph *sigVsPres = new TGraph(nfiles, pressures, means);
  Int_t graphPadNum = nfiles+1;
  cout << "graph pad num = " << graphPadNum << endl;
  c_pres->cd(graphPadNum);
  sigVsPres->Draw("APL");
  sigVsPres->GetXaxis()->SetTitle("Pressure (torr on convectron)");
  sigVsPres->GetYaxis()->SetTitle("MCA peak voltage (mV)");
  c_pres->Update();

  c_pres->Print("heightVsPressure.pdf");

  TFile outfile("zcoord_analysis.root", "RECREATE");
  sigVsPres->Write();
  outfile.Close();

}

void mca(TString fname) {

  //c_scope = new TCanvas("scopemca","mca",400,400); 

  d = new DmtpcDataset;
  d->openRootFile(fname);

  nEntries = d->tree()->GetEntries();
  cout << "d->tree()->GetEntries() = " << d->tree()->GetEntries() << endl;

  // FIX hist range...
  TH1F *mca1 = new TH1F("mca1", "Scope MCA", 100, -20., 20.); // millivolts

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

  mca1->Fit("gaus", "", "", vfitMin[mcaNumber], vfitMax[mcaNumber]);
  TF1 *fit = mca1->GetFunction("gaus");
  Float_t mean = fit->GetParameter(1);
  means[mcaNumber] = mean;
  cout << "mean = " << mean << endl;
  cout << "there were " << nTrigTotal << " triggers" << endl;

  //TFile outfile("mcaHist.root", "RECREATE");
  //mca1->Write();
  //outfile.Close();
}

