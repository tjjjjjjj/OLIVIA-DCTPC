//#include "TTimeStamp.h"
#include "TDatime.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "../MaxCamChannel.hh"
using namespace std;

//int logPressure() {
//  cout << "hi" << endl;
//  return 0;
//}

int logPressure(int nsleep=500) {
  // inputs:
  //    nsleep  --  milliseconds between reads
  //
  TFile *file = new TFile("pressure.root", "RECREATE", "example root file");
  TTree *tree = new TTree("tree", "A ROOT Tree");

  //TTimeStamp *time = new TTimeStamp;
  TDatime *time = new TDatime;
  Float_t prsCDG, prsBPG;
  
  //tree->Branch("time", "TTimeStamp", &time);
  tree->Branch("time", "TDatime", &time);
  tree->Branch("prsCDG", &prsCDG, "prsCDG/F");
  tree->Branch("prsBPG", &prsBPG, "prsBPG/F");

  char *inficonPort="/dev/ttyS1";
  Int_t cdgChannel = 2;
  Int_t bpgChannel = 1;

  MaxCamChannel cdg("cdg", "cdg");
  MaxCamChannel bpg("bpg", "bpg");

  //Int_t nmax = 2;
  //for (Int_t ii=0; ii<nmax; ii++) {
  while(1) {
    time->Set();
    prsCDG = cdg.readInficonController(inficonPort, cdgChannel);
    prsBPG = bpg.readInficonController(inficonPort, bpgChannel);
    cout << "time, prsCDG, BPG = " << time->Get() 
	 << ", " << prsCDG << ", " << prsBPG << endl;
    tree->Fill();
    // force a write to disk
    tree->AutoSave("SaveSelf");
    // Flush doesn't work...?
    //tree->FlushBaskets();  // write to disk...

    gSystem->Sleep(nsleep);
  }

  file->Write();
  //tree->Print();
  return 0;
}
