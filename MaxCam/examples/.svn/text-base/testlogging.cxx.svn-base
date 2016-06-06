// this code was used to get the logging off of the ground and to
// solve a problem with writing/reading objects in TClonesArrays to
// disk.
// See DetectorLogger.cxx for the real code though...
//

#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "../DmtpcLoggableParam.hh"

const int MAX_ENTRIES = 100;
const int NLOOP = 15;
const int SPLIT = 0;

int debug = 0;  // can be 0, 1, 2, 3, ... for various levels of debug

TClonesArray *pres_tca_ptr = new TClonesArray("DmtpcLoggableParam", MAX_ENTRIES);
TClonesArray &pres_tca = *pres_tca_ptr;
TClonesArray *temp = new TClonesArray("DmtpcLoggableParam", MAX_ENTRIES);

DmtpcLoggableParam presParam("pressure0", "cf4 pressure");
DmtpcLoggableParam tempParam("temperature0", "cf4 temperature");

TFile *file;
TTree *tree;
Int_t nFills = 0;
Int_t nPres = 0;
Int_t nTemp = 0;
Int_t pressure_interval = 2;
Int_t temperature_interval = 3;
Int_t filltree_interval = 5;
TString filename("logtest.root");

void read_pressure(int ii);
void read_temperature(int ii);
void reset();
void reset_pressure();
void reset_temperature();
void fill_tree();

void testlogging() {
  cout << "filename = " << filename << endl;
  file = new TFile(filename, "RECREATE");
  tree = new TTree("env", "env");
  tree->Branch("pressure", "TClonesArray", &pres_tca_ptr, 32000, SPLIT);
  tree->Branch("temperature", "TClonesArray", &temp, 32000, SPLIT);

  // not sure exactly what this means...
  //pres_tca_ptr->BypassStreamer();

  reset();
  for (int ii=0; ii<NLOOP; ii++) {
    cout << "ii = " << ii << endl;
    if (ii%pressure_interval == 0) {
      read_pressure(ii);
    }
    //if (ii%temperature_interval == 0) {
    //  read_temperature(ii);
    //}
    if (ii%filltree_interval == 0) {
      fill_tree();
      cout << "done with fill number " << nFills << endl;
      cout << "" << endl;
    }
  }
  // then, at the end, fill the tree to capture whatever data has not
  // yet been logged...
  fill_tree();


  //file = tree->GetCurrentFile();
  if (debug > 3) cout << "here3" << endl;
  tree->Print();
  tree->Write();
  //file->Write();
  if (debug > 3) cout << "here4" << endl;
  file->Close();
  if (debug > 3) cout << "here5" << endl;

  nFills = 0;
}

void fill_tree() {
  cout << "fill_tree()" << endl;
  tree->Fill();
  nFills++;
  reset();
}

void reset() {
  reset_pressure();
  reset_temperature();
}

void reset_pressure() {
  nPres = 0;
  // clear the TClonesArray in preparation for the next fill
  pres_tca.Clear();
  //pres_tca_ptr->Delete();
}
void reset_temperature() {
  nTemp = 0;
  // clear the TClonesArray in preparation for the next fill
  temp->Delete();
}

void read_pressure(int ii) {
  cout << "read_pressure()" << endl;
  presParam.SetValue((ii+1)*3.0);
  presParam.SetTimeStamp(ii);
  new ( pres_tca[nPres] ) DmtpcLoggableParam(presParam);
  cout << "echo:  pres->At("<<nPres<<") holds " 
       << ((DmtpcLoggableParam*)pres_tca_ptr->At(nPres))->GetTimeStamp().GetTime() << ", "
       << ((DmtpcLoggableParam*)pres_tca_ptr->At(nPres))->GetValue() << endl;
 nPres++;
}

void read_temperature(int ii) {
  cout << "read_temperature()" << endl;
  tempParam.SetValue((ii+1)*7.0);
  tempParam.SetTimeStamp(ii);
  new ( (*temp)[nTemp] ) DmtpcLoggableParam(tempParam);
  cout << "echo:  temp->At("<<nTemp<<") holds " 
       << ((DmtpcLoggableParam*)temp->At(nTemp))->GetTimeStamp().GetTime() << ", "
       << ((DmtpcLoggableParam*)temp->At(nTemp))->GetValue() << endl;
  nTemp++;
}

TTree *readTree() {
  cout << "filename = " << filename << endl;
  TFile *fin = new TFile(filename);
  TTree *tin = (TTree*)fin->Get("env");
  return tin;
}

void read() {
  cout << "filename = " << filename << endl;
  TFile *fin = new TFile(filename);
  TTree *tin = (TTree*)fin->Get("env");
  
  cout << "tin->GetEntries() = " << tin->GetEntries() << endl;

  TClonesArray *in_pres_tca = new TClonesArray("DmtpcLoggableParam", MAX_ENTRIES);
  tin->GetBranch("pressure")->SetAutoDelete(kFALSE);  // ?????what does this do?????
  tin->SetBranchAddress("pressure", &in_pres_tca);

  TClonesArray *temp_tca = new TClonesArray("DmtpcLoggableParam", MAX_ENTRIES);
  tin->SetBranchAddress("temperature", &temp_tca);

  for (int ii=0; ii<tin->GetEntries(); ii++) {
    in_pres_tca->Clear();
    cout << "Event " << ii << endl;
    tin->GetEntry(ii);
    cout << "----- in_pres_tca ----- " << endl;
    //cout << "in_pres_tca->Print()" << endl;
    //in_pres_tca->Print();
    for (int jj=0; jj<in_pres_tca->GetEntries(); jj++) {
      DmtpcLoggableParam *param = (DmtpcLoggableParam*)in_pres_tca->At(jj);
      cout << "  jj, name, time, pres = " << jj << ", " 
	   << param->GetParamName() << ", " 
	   << param->GetTimeStamp().GetTime() << ", " 
	   << param->GetValue() << endl;
      //<< ((DmtpcLoggableParam*)in_pres_tca->At(jj))->GetValue() << endl;
    }
    cout << "" << endl;
    cout << "----- temp_tca ----- " << endl;
    //cout << "temp_tca->Print()" << endl;
    //temp_tca->Print();
    for (int jj=0; jj<temp_tca->GetEntries(); jj++) {
      cout << "  jj, temp = " << jj << ", " 
    	   << ((DmtpcLoggableParam*)temp_tca->At(jj))->GetValue() << endl;
    }
    cout << "" << endl;
    cout << "" << endl;
  }

  fin->Close();
}
