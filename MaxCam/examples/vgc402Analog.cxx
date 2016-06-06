// Add automatic filename generation  YYMMDD_log.root (perhaps?)
// one file per day
// changeover happens at midnight UTC (perhaps?)
// times should all be stored as UTC (*not* local time!)
// In while loop, need a test to see if file should be closed and new file opened...
// So we need an init() section that is called when prog is run first time
// and then is called once per day in the loop

#include <fstream>
//#include <string>
#include <iostream>
//#include "TDatime.h"
#include "TTimeStamp.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "../MaxCamChannel.hh"

using namespace std;

void vgc402Analog(char *outputfile="chamberPumpdown.dat") {

  //char *outputfile = "chamberPumpdown.dat";
  ofstream out;
  out.open(outputfile);
  TString hdr = "# yyyy-mm-dd hh:mm:ss secSince1995 BPG400_pressure_volts";
  out  << hdr << endl;
  cout << hdr << endl;
  
  TFile *file = new TFile("chamberPumpdown.root", "RECREATE", "example root file");
  TTree *tree = new TTree("tree", "A ROOT Tree");
  
  Int_t sleepTimeMilliseconds = 5000;
  TTimeStamp *time = new TTimeStamp;
  Float_t prsBPG;
  UInt_t timeUInt;

  tree->Branch("prsBPG", &prsBPG, "prsBPG/F");
  tree->Branch("time", &timeUInt, "timeUInt/i");

  while(1) {
    time->Set();
    MaxCamChannel *bpg = new MaxCamChannel("chanA", "BPG", 0);
    prsBPG = bpg->ni6229ReadADC();
    timeUInt = UInt_t(time->AsDouble()); // TTimeStamp seconds since 1970 Jan 1 

    cout << time->AsString("s") << " " << timeUInt << "  " << prsBPG << endl;
    out  << time->AsString("s") << " " << timeUInt << "  " << prsBPG << endl;

    tree->Fill();
    tree->AutoSave("SaveSelf");

    gSystem->Sleep(sleepTimeMilliseconds);
  }

  file->Write();

}
