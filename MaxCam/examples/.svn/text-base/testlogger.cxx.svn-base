#include <iostream>
#include <vector>
#include <unistd.h> // usleep()

#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
//#include "TNamed.h"
#include "TTimeStamp.h"

#include "../DmtpcLogEntry.hh"
#include "../DmtpcLoggableParam.hh"

using namespace std;

#define SPLIT 0

void testlogger() {
  TString logfilename_full("junk.root");
  TFile *logfile = new TFile(logfilename_full, "RECREATE", "example log file");
  TTree *logtree = new TTree("environment", "environment");

  DmtpcLogEntry *presCF4 = new DmtpcLogEntry("presCF4", "CDG gauge");
  logtree->Branch("presCF4", "DmtpcLogEntry", &presCF4, 32000, SPLIT);
  
  DmtpcLoggableParam *paramCF4 = new DmtpcLoggableParam("presCF4", "cdg gauge reading");
  TTimeStamp *timeCF4 = new TTimeStamp();

  for (Int_t ii=0; ii<8; ii++) {
    cout << "ii = " << ii << endl;
    if (ii%2 == 0) {
      cout << "timeCF4->Set()" << endl;
      timeCF4->Set();
      cout << "paramCF4->SetValue()" << endl;
      paramCF4->SetValue(ii+1);
      cout << "paramCF4->SetTimeStamp()" << endl;
      paramCF4->SetTimeStamp(timeCF4);
      cout << "presCF4->Append()" << endl;
      presCF4->Append(paramCF4);
    }

    if (ii%3 == 0) {
      // fill tree
      cout << "pre-fill:  presCF4 holds N elements:  N=" << presCF4->GetNSamples() << endl;
      cout << "------------------- fill --------------------- " << endl;
      logtree->Fill();
      cout << "post-fill nelements = " << presCF4->GetNSamples() << endl;
      cout << "presCF4->Reset()" << endl;
      presCF4->Reset();
      cout << "presCF4->Reset() done" << endl;
    }
    usleep(1000000);
  }

  logtree->Fill();
  presCF4->Reset();
  
  logtree->Print();
  logtree->Write();
  logfile->Close();

}


void drawloop() {
  TFile f("junk.root");
  TTree *t = (TTree*)f.Get("environment");

  TGraph *gPresCF4 = new TGraph();

  DmtpcLogEntry *dle = new DmtpcLogEntry();
  cout << "setbranchaddress" << endl;
  t->SetBranchAddress("presCF4", &dle);
  Int_t nPts = 0;
  for (int ii=0; ii<3; ii++) {
    t->GetEvent(ii);
    cout << dle->GetNSamples() << endl;
    dle->print();
    for (Int_t jj=0; jj<dle->GetNSamples(); jj++) {
      gPresCF4->SetPoint(nPts, 
			 dle->GetTimeStamp(jj)->GetTime(), 
			 dle->GetValue(jj));
      nPts++;
    }
  }

  gPresCF4->Draw("alp");
}

void draw() {
  TFile f("junk.root");
  TTree *t = (TTree*)f.Get("environment");
  //t->Draw("presCF4->GetNSamples()");
  //t->Draw("presCF4->GetValue()");  // doesn't work
  t->Draw("presCF4->GetSamples()");  // doesn't work
  
}

void view() {
  TFile f("junk.root");
  TTree *t = (TTree*)f.Get("environment");
  t->Scan("((TClonesArray*)(presCF4@.GetSamples())).GetValue()");
}

void read2() {
  TFile f("junk.root");
  TTree *t = (TTree*)f.Get("environment");

  DmtpcLogEntry *dle = new DmtpcLogEntry();
  TClonesArray *tca = new TClonesArray("DmtpcLoggableParam");
  cout << "setbranchaddress" << endl;
  t->SetBranchAddress("presCF4", &dle);
  for (int ii=0; ii<3; ii++) {
    t->GetEvent(ii);
    cout << dle->GetNSamples() << endl;
    tca = dle->GetSamples();
    cout << "nentries = " << tca->GetEntries() << endl;
    dle->print();
  }
}

void read() {
  TFile f("junk.root");
  TTree *t = (TTree*)f.Get("environment");

  DmtpcLogEntry *dle = new DmtpcLogEntry();
  cout << "setbranchaddress" << endl;
  t->SetBranchAddress("presCF4", &dle);
  Int_t nentries = t->GetEntries();
  cout << "nentries = " << nentries << endl;
  for (Int_t ii=0; ii<nentries; ii++) {
    cout << "ii = " << ii << endl;
    cout << "getevent" << endl;
    t->GetEvent(ii);
    cout << "after getevent" << endl;
    //dle->print();
    for (Int_t jj=0; jj<dle->GetNSamples(); jj++) {
      cout << "jj = " << jj << endl;
      cout << "time, val = " << dle->GetTimeStamp(jj)->AsString("s") << ", " << dle->GetValue(jj) << endl;
    }
    cout << endl;
  }


}

//void testme() {
//  
//  TFile f("junk.root");
//  TTree *t = (TTree*)f.Get("environment");
//  DmtpcLogEntry *dle = new DmtpcLogEntry();
//  t->SetBranchAddress("presCF4", &dle);
//  t->GetEvent(2);
//
//  vector<Float_t> myvals;
//  cout << myvals.size() << endl;
//  dle->GetValues(myvals);
//  cout << myvals.size() << endl;
//  for (uint ii=0; ii<myvals.size(); ii++) {
//    cout << myvals[ii] << " ";
//  }
//  cout << endl;
//}

//void unused() {
//
//  TClonesArray *tca_ptr = new TClonesArray("TNamed");
//  TNamed mynamed("naem", "desc");
//  new ( (*tca_ptr)[0] ) TNamed(mynamed);
//  cout << tca_ptr->At(0)->GetName() << endl;
//
//  //TClonesArray *tca2_ptr = new TClonesArray("TNamed");
//  TClonesArray *tca2_ptr;
//  tca2_ptr = (TClonesArray*)tca_ptr->Clone();
//
//  cout << tca2_ptr->At(0)->GetName() << endl;
//
//  ((TNamed*)tca2_ptr->At(0))->SetName("changed");
//  
//  cout << tca_ptr->At(0)->GetName() << endl;
//  cout << tca2_ptr->At(0)->GetName() << endl;
//}
