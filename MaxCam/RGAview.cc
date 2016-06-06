#include "TObjArray.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjString.h"
#include "TTimeStamp.h"
#include "TGraph.h"
#include "TH2F.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "RGAview.hh"
#include "TVectorD.h"

using namespace std;

ClassImp(RGAview)

RGAview::RGAview(char RGAdata[50], char filename[50]) {
  strcpy(datafile, RGAdata);
  strcpy(rootfile, filename);
  scan_number = 1;
  mass_number = 1;
  nmasses = 0;
}

RGAview::RGAview(const RGAview & object) {
  strcpy(datafile, object.datafile);
  strcpy(rootfile, object.rootfile);
  scan_number = 1;
  mass_number = 1;
  nmasses = 0;
}

RGAview::~RGAview() {
}

int RGAview::readRGA() {

  //First, we read out the header, which contains the info about the masses
  Double_t masses[50] = {0};
  nmasses = 50;
 
  readheader(masses);
  
  //We then read the pressure scans and time stamps into the tree
  readdata(masses);
  return 1;
}

void RGAview::readheader(Double_t * massarray) {
  ifstream infile;
  string file = datafile;
  infile.open(file.c_str());
  if(!infile) {
    cout << "Cannot open file!" << endl;
  }
  //first, we clean up the header
  string header;
  //Skip the blank line at the beginning
  getline(infile, header);
  getline(infile, header);
  TString Theader = (TString)header;
  Theader.ReplaceAll("\"", "");
  Theader.ReplaceAll("Mass ", "");  
  Theader.ReplaceAll("Scan\t", "");  
  Theader.ReplaceAll("Time\t", "");
  Theader.ReplaceAll(".0", "");
  Theader.ReplaceAll("\tSum Scanned Masses\t", "");

  //Then, we read the header into massarray
  TObjArray *mass_obj_array = Theader.Tokenize("\t");
  Int_t ii;
  Int_t nmasses = mass_obj_array->GetEntries();
  for(ii = 0; ii < nmasses; ii++) {
    //convert form TObjString->TString->float...whew!
    massarray[ii] = ((TObjString*)mass_obj_array->At(ii))->GetString().Atof();
  }
  infile.close();
}

void RGAview::readdata(Double_t * massarray) {
  ifstream infile;
  string file = datafile;
  infile.open(file.c_str());
  if(!infile) {
    cout << "Cannot open file!" << endl;
  }
  //We convert each scan into a TGraph and place it inside of "RGAtree"
  string oneFullScan;
  TString scan;
  TString scan_time_stamp;
    
  //We begin by creating the file and tree that will ultimately store everything
  string rfile = rootfile;
  TFile *t_file = new TFile(rfile.c_str(), "RECREATE");
  TTree *RGA_tree = new TTree("RGA_data", "RGA_data");

  Double_t pressures[50] ={0};
  TGraph *graph_scan = new TGraph(nmasses, massarray, pressures);
  RGA_tree->Branch("PressureScans", "TGraph", &graph_scan);

  TDatime *time = new TDatime();
  RGA_tree->Branch("TimeStamp", "TDatime", &time);
  Double_t totalpressure; 
  RGA_tree->Branch("TotalPressure", &totalpressure, "totalpressure/D"); 

  //skip the header line and the blank line
  getline(infile, oneFullScan);
  getline(infile, oneFullScan);

  // loop over each line in file untill the end is reached
  getline(infile, oneFullScan);
  while (!infile.eof()) {
    scan = TString(oneFullScan);
    scan.ReplaceAll("\"","");
    TObjArray *tokens = scan.Tokenize("\t");

    //extract the time stamp
    scan_time_stamp = ((TObjString*)tokens->At(0))->GetString();
    convert_to_TDatime(time, scan_time_stamp);

    //Extract the pressure data
    for(Int_t ii=2; ii<52; ii++) {
      pressures[ii-2] = ((TObjString*)tokens->At(ii))->GetString().Atof();
      graph_scan->SetPoint(ii-2, massarray[ii-2], pressures[ii-2]);
    } 
    
    totalpressure = 0;
    for(Int_t i=0; i<50; i++) {
      totalpressure += pressures[i];
    }
   
    //Write to the tree
    RGA_tree->Fill();
    delete tokens;
    getline(infile, oneFullScan);
  }
  infile.close();
  //Finally, we write the tree to t_file and close the TFile
  RGA_tree->Write();
  t_file->Close();
}

void RGAview::convert_to_TDatime(TDatime *ttime, TString time_str) {
  time_str.ReplaceAll("-"," ");
  time_str.ReplaceAll(":"," ");
  time_str.ReplaceAll("."," ");
  time_str.ReplaceAll("/"," ");
  TObjArray *time_obj = time_str.Tokenize(" ");
  Int_t year = ((TObjString*)time_obj->At(0))->GetString().Atoi();
  Int_t month = ((TObjString*)time_obj->At(1))->GetString().Atoi();
  Int_t day = ((TObjString*)time_obj->At(2))->GetString().Atoi();
  Int_t hours = ((TObjString*)time_obj->At(3))->GetString().Atoi();
  Int_t minutes = ((TObjString*)time_obj->At(4))->GetString().Atoi();
  Int_t seconds = ((TObjString*)time_obj->At(5))->GetString().Atoi();
  *ttime = TDatime(year, month, day, hours, minutes, seconds);
}

void RGAview::view_one_scan(Int_t n) {
  scan_number = n;
  string file = rootfile;
  TFile *in = new TFile(file.c_str());
  TTree *alldata = (TTree*)in->Get("RGA_data");
  TGraph *onescan = new TGraph();
  alldata->SetBranchAddress("PressureScans", &onescan);
  alldata->GetEvent(n);
  string scan;
  stringstream out;
  out << n;
  scan = out.str();
  string title = "Scan " + scan;;
  onescan->SetTitle(title.c_str());
  onescan->Draw("alp");
  onescan->GetHistogram()->GetXaxis()->SetTitle("Mass (AMU)");
  onescan->GetYaxis()->SetTitle("Pressure (torr)");
  onescan->GetXaxis()->CenterTitle();
  onescan->GetYaxis()->CenterTitle();
}

void RGAview::next_scan() {
  scan_number++;
  view_one_scan(scan_number);
}

void RGAview::view_one_mass(int n) {
  //First we extract the data from the saved tree
  mass_number = n;
  string file = rootfile;
  TFile *in = new TFile(file.c_str());
  TTree *alldata = (TTree*)in->Get("RGA_data");
  TGraph *onescan = new TGraph();
  alldata->SetBranchAddress("PressureScans", &onescan);
  TDatime *onetime = new TDatime();
  alldata->SetBranchAddress("TimeStamp", &onetime);

  //Next, we create a TGraph that will contain the data for mass n
  TGraph *onemass = new TGraph();
  
  //Then, we copy the desired data from "alldata" into "onemass"
  Int_t num_entries = alldata->GetEntries();
  for(int i=0; i<num_entries; i++) {
    alldata->GetEvent(i);
    Double_t *array = onescan->GetY();
    Double_t pressure = array[n-1];
    onemass->SetPoint(i, onetime->Convert(), pressure); 
  }

  //Finally, we plot onemass
  string mass;
  stringstream out;
  out << n;
  mass = out.str();
  string title = "Mass " + mass;
  onemass->SetMarkerStyle(21);
  onemass->SetMarkerSize(1);
  onemass->GetXaxis()->SetTimeDisplay(1);
  onemass->GetXaxis()->SetTimeOffset(0, "gmt");
  onemass->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M:%S");
  onemass->GetXaxis()->SetTitle("Time");
  onemass->GetYaxis()->SetTitle("Pressure (torr)");
  onemass->GetXaxis()->CenterTitle();
  onemass->GetYaxis()->CenterTitle();
  onemass->GetXaxis()->SetLabelSize(0.0125);
  onemass->SetTitle(title.c_str());
  onemass->Draw("ap");
}

void RGAview::next_mass() {
  mass_number++;
  view_one_mass(mass_number);
}

void RGAview::view_total_pressure() {
  //First we extract the data from the saved tree
  string file = rootfile;
  TFile *in = new TFile(file.c_str());
  TTree *alldata = (TTree*)in->Get("RGA_data");
  TDatime *onetime = new TDatime();
  alldata->SetBranchAddress("TimeStamp", &onetime);
  Double_t total_pressure;
  alldata->SetBranchAddress("TotalPressure", &total_pressure);

  //Next, we create a TGraph that will contain the data for mass n
  TGraph *pressure = new TGraph();
  
  //Then, we copy the desired data from "alldata" into "pressure"
  Int_t num_entries = alldata->GetEntries();
  for(int i=0; i<num_entries; i++) {
    alldata->GetEvent(i);
    pressure->SetPoint(i, onetime->Convert(), total_pressure); 
  }

  //Finally, we plot onemass
  string title = "Total Pressure vs. Time";
  pressure->SetMarkerStyle(21);
  pressure->SetMarkerSize(1);
  pressure->GetXaxis()->SetTimeDisplay(1);
  pressure->GetXaxis()->SetTimeOffset(0, "gmt");
  pressure->GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M:%S");
  pressure->GetXaxis()->SetTitle("Time");
  pressure->GetYaxis()->SetTitle("Pressure (torr)");
  pressure->GetXaxis()->CenterTitle();
  pressure->GetYaxis()->CenterTitle();
  pressure->GetXaxis()->SetLabelSize(0.0125);
  pressure->SetTitle(title.c_str());
  pressure->Draw("apl");
}
