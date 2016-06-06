#include "TObjArray.h"
#include "TH1F.h"
#include "TString.h"
//
#include <string>
#include <iostream>
#include <fstream>
//
using namespace std;

// read and display RGA scans
// assumes that you have dummy RGA data (see rgadummy.dat)
//
// > root -l
// root [0] .L readRGA.cxx+
// root [1] go()
// root [2] view(0)  // view the first scan
// root [3] view(1)  // view the second scan

// an array of objects to hold all of the RGA scans
// A TObjArray can hold different objects.  In this
// case we are only putting TH1F objects into the array
TObjArray* allScans = new TObjArray();

int go() {
  // this function will read and store in memory the RGA scans
  // eventually, you want to read the ASCII RGA data and 
  // write to a ROOT data tree/file
  // Each entry in the tree will contain a TH1F of an RGA scan and a timestamp

  // open up a file that contains dummy RGA data
  TString filename = "rgadummy.dat";
  ifstream infile;
  infile.open(filename);
  if (!infile) {
    cout << "cannot open file." << endl;
    return 1;
  }

  // bin range definition for the TH1F to store a single RGA scan
  Int_t massMin = 1;
  Int_t massMax = 100;
  Int_t nmasses = 100;

  string oneFullScan;
  TString scan;

  Int_t scanNumber = 0;

  // loop over each line in the file until you reach the end of file
  while (! infile.eof() ) {
    // read a full line of the data file
    // oneFullScan should hold "1 2 3 4 5 6 7 ... 100"
    getline(infile, oneFullScan);
    cout << oneFullScan << endl;

    // ROOT's TString has a nice, clean way to break a string
    // into tokens based on a delimiter
    TString scan(oneFullScan);
    TObjArray *tokens = scan.Tokenize(" ");
    // tokens now holds ["1", "2", "3", ..., "100"]
    // but each entry is a TObjString, not just a TString

    Int_t ntokens = tokens->GetEntries();
    TString histName = "h_scan";
    histName += scanNumber;
    TH1F *h_scan = new TH1F(histName, histName, nmasses, massMin, massMax+1);
    Float_t partialPressure;
    // now, loop through each token in the TObjArray, 
    // get the TObjString, then get the TString, then convert to a float
    // finally, store that value into the histogram
    for (Int_t ii=1; ii<=ntokens; ii++) {
      partialPressure = ((TObjString*)tokens->At(ii-1))->GetString().Atof();
      h_scan->SetBinContent(ii, partialPressure);
    }
    // add this scan into the array that holds all scans
    allScans->Add(h_scan);
    delete tokens;
    scanNumber++;
  }

  return 0;
}

void view(Int_t nscan) {
  // after scans have been read from ASCII file into memory
  // you can view the scan with this function
  // Argument is the scan number to view (starting at zero)
  ((TH1F*)allScans->At(nscan))->Draw();
}
