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
int iev;
TCanvas *c_ccd=0, *c_scope=0;
TString what="ccd scope";

const int maxWF=8;

void init(TString fname) {
    if (what.Contains("ccd")) { c_ccd= new TCanvas("c_ccd","", 600, 400); c_ccd->Divide(2,1); }
    if (what.Contains("scope")) { c_scope= new TCanvas("c_scope","",400,800); c_scope->Divide(2,maxWF/2); }
    d= new DmtpcDataset;
    d->openRootFile(fname);
    for (int iccd=0; iccd<2; iccd++) {
        bias[iccd]=d->getBiasFrame(iccd+1);
    }
    iev=0;
}


void next() {
  d->getEvent(iev++);
  if (what.Contains("ccd")) {
      for (int iccd=0; iccd<d->event()->ccdData()->GetEntries(); iccd++) {
          c_ccd->cd(iccd+1);
          image[iccd]=d->event()->ccdData(iccd);
          if (!image[iccd]) continue;
          image[iccd]->Add(bias[iccd],-1);
          image[iccd]->SetMinimum(-50);
          image[iccd]->SetMaximum(250);
          image[iccd]->DrawCopy("colz");
      }
  }
  if (what.Contains("scope")) {
      int nwf=d->event()->scopeData()->GetEntries()/2;
      int iwf=0;
      for (int i=0; i<nwf; i++) {
          if (2*iwf>maxWF) break;
          c_scope->cd(2*iwf+1); d->event()->scopeData(i)->Draw();
          c_scope->cd(2*iwf+2); d->event()->scopeData(i+nwf)->Draw();
          MaxCamWaveformTools wfA(d->event()->scopeData(i)); wfA.print();
          MaxCamWaveformTools wfB(d->event()->scopeData(i+nwf)); wfB.print();
          iwf++;
      }
      cout << "total wf="<<nwf << " plotted="<<iwf<<endl;
  }
}

void loop(int n) {
  for (int i=0; i<n; i++) { 
    cout << "event="<<iev << endl;
    next(); 
  }
}
