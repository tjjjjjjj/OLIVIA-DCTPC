//   root[0] .L process.cxx++
//

//#include <exception>

// confirmed
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TArrayF.h"
#include "TSystem.h"
#include <string>
#include "TObjArray.h"
#include "TObjString.h"

// leftovers
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "../DmtpcDataset.hh"
#include "../MaxCamWaveformTools.hh"
#include "../ScopeDataInfo.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
//using std::cout;
//using std::endl;
using namespace std;

DmtpcDataset *d;
TH2F *image[2], *bias[2];
int iev, nEntries;
bool firstPass;
// # of channels on the scope.  
// this matters because the wfs from all scope channels 
// are returned in a 1d array
Int_t nChan = 2;  

bool fileIsValid(TString filename);
void processFile(TString infile);
void processFileRange(TString fileroot, Int_t min, Int_t max);
Float_t calcMCAMean(Int_t runnum);
void getMCAMean(TString filename);
//TString fileOfRun(Int_t runnum, TString fileroot="zcoord_VS_B_");
TString fileOfRun(Int_t runnum, TString fileroot="/data/zcoord/zcoord_VS_E_");
void viewMCA(Int_t runnum);
TH1F* makeMCA(Int_t runnum);
TH1F* makeMCA(TString filename);
void plotMeans(Int_t sleep_ms=0);
TF1* fitToMCA(TH1F* mca);
void processListOfFiles(TString filename="overnightFiles.dat");
TObjArray* getFilenames(TString listname, Bool_t verbose=false);
void plotGainVsTime(TString filename="overnightFiles.dat", Int_t sleep_ms=0);

const Int_t ndata = 13;
const float motorHeights[ndata] = {0., 5, 10, 15, 20, 25, 30, 25, 20, 15, 10, 5, 0};
TArrayF heights = TArrayF(ndata, motorHeights);

Int_t pAmpl = 0;
Int_t pMean = 1; // the second parameter in gaus
Int_t pStdDev = 2;

TCanvas *c_gr  = new TCanvas("gr","gr");
TCanvas *c_mca = new TCanvas("mca", "mca");
TCanvas *c_means = new TCanvas("means", "means");

void testgraph() {
  const Int_t n = 10;
  Double_t x[n]  = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
  Double_t y[n]  = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
  Double_t ex[n] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
  Double_t ey[n] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
  TGraphErrors *gr = new TGraphErrors(n,x,y,ex,ey);
  gr->SetTitle("TGraphErrors Example");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("ALP");
}

void plotGainVsTime(TString filename, Int_t sleep_ms) {
  TObjArray* filelist = getFilenames(filename);
  TIter next(filelist);
  Int_t nfiles = filelist->GetEntries();
  TArrayF means(nfiles);
  TArrayF meanErrs(nfiles);
  TArrayF xvals(nfiles);
  c_means->cd();
  c_means->Divide(4,4);
  Int_t ican = 1;
  while (TObjString* rootFile=(TObjString*)next()) {
    c_means->cd(ican);
    TString filename = rootFile->GetString();
    filename.ReplaceAll(".root", "_red.root");
    cout << "filename = " << filename << endl;
    TH1F* mca = makeMCA(filename);
    TF1 *f1   = fitToMCA(mca);
    means.AddAt(-1000.*f1->GetParameter(pMean), ican-1);
    //meanErrs.AddAt(1000*mca->GetRMS(), ican-1);
    meanErrs.AddAt(TMath::Abs(1000*f1->GetParameter(pStdDev)), ican-1);
    xvals.AddAt(ican*1.0, ican-1);
    ican++;
    c_means->Update();
    gSystem->Sleep(sleep_ms);
    delete mca;
    delete f1;
  }

  for (Int_t ii=0; ii<nfiles; ii++) {
    cout << xvals.At(ii) << "  " << means.At(ii) << endl;
  }

  c_gr->cd();
  TGraphErrors *gg = new TGraphErrors(nfiles, xvals.GetArray(), means.GetArray(), NULL, meanErrs.GetArray());
  gg->Draw("AL*");
  gg->GetHistogram()->SetXTitle("Time (arb units)");
  gg->GetHistogram()->SetYTitle("Mean pulse height [mV]");

  TF1* fitFunc = new TF1("fitFunc", "[0]+[1]*exp(-(x-[2])/[3])", 0, 80);
  //TF1* fitFunc = new TF1("fitFunc", "[0]+expo(1)", 0, 80);
  fitFunc->SetParameter(0, 10.0);
  fitFunc->SetParameter(1, 50.0);
  fitFunc->SetParameter(2,  1.0);
  fitFunc->SetParameter(3, 20.0);
  gg->Fit("fitFunc");
  c_gr->Update();

  TString outfile_graph("graph.root");
  TFile graphFile(outfile_graph, "RECREATE");
  gg->Write();
  graphFile.Close();

  // clean up
  delete filelist;
  //delete gg;
}

void plotMeans(Int_t sleep_ms) {
  // should pass fileroot and use
  // TString fileOfRun(Int_t runnum, TString fileroot) instead of hardcoding..
  TArrayF means(ndata);
  TArrayF meanErrs(ndata);
  c_means->cd();
  c_means->Divide(4,4);
  for (Int_t ii=1; ii<=ndata; ii++) {
    c_means->cd(ii);
    TH1F *mca = makeMCA(ii);
    TF1 *f1   = fitToMCA(mca);
    means.AddAt(-1000.*f1->GetParameter(pMean), ii-1);
    //meanErrs.AddAt(1000*f1->GetParError(pMean), ii-1);
    meanErrs.AddAt(TMath::Abs(1000*f1->GetParameter(pStdDev)), ii-1);
    //meanErrs.AddAt(1000*mca->GetRMS(), ii-1);
    gSystem->Sleep(sleep_ms);
  }
  c_means->Update();

  c_gr->cd();
  TGraphErrors *gg = new TGraphErrors(ndata, heights.GetArray(), means.GetArray(), NULL, meanErrs.GetArray());
  gg->Draw("AL*");
  gg->GetHistogram()->SetXTitle("Motor position (raw) [mm]");
  gg->GetHistogram()->SetYTitle("Mean pulse height [mV]");
  c_gr->Update();
}

void viewMCA(Int_t runnum) {
  TString filename = fileOfRun(runnum);
  filename.ReplaceAll(".root", "_red.root");
  cout << "filename = " << filename << endl;

  TFile *ff = new TFile(filename);
  TTree *tt = (TTree*)ff->Get("outtree");
  c_mca->cd();
  tt->Draw("pulseHeight");
  c_mca->Update();
}

void getMCAMean(TString filename) {
  // compute pulse heights and store in a new root file
  processFile(filename);
  TString redfile = TString(filename);
  redfile.ReplaceAll(".root", "_red.root");

  TH1F *mca = makeMCA(redfile);
  TF1 *f1   = fitToMCA(mca);
  Float_t mean = -1000.*f1->GetParameter(pMean);
  Float_t std  = 1000*f1->GetParameter(pStdDev);
  cout << "mean, stdev = " << mean << ", " << std << endl;
  
}

TH1F* makeMCA(TString filename) {
  
  TFile *ff = new TFile(filename);
  TTree *tt = (TTree*)ff->Get("outtree");
  tt->Draw("pulseHeight");
  TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
  htemp->SetXTitle("Pulse Height [V]");
  return htemp;
}

TH1F* makeMCA(Int_t runnum) {
  TString filename = fileOfRun(runnum);
  filename.ReplaceAll(".root", "_red.root");
  cout << "filename = " << filename << endl;
  return makeMCA(filename);
}

TF1* fitToMCA(TH1F* mca) {
  //TF1 *fitFunc = new TF1("fitFunc", "gaus(0)+[3]", -0.1, 0.0);
  TF1 *fitFunc = new TF1("fitFunc", "gaus(0)", -0.1, 0.0);
  cout << "mca->GetMaximum() = " << mca->GetMaximum() << endl;
  cout << "mca->GetMean() = " << mca->GetMean() << endl;
  cout << "mca->GetBinContent(mca->GetMaximumBin()) = " << mca->GetBinContent(mca->GetMaximumBin()) << endl;
  //Float_t mcaMean = mca->GetMean();
  Float_t mcaPeak = mca->GetBinContent(mca->GetMaximumBin());
  Float_t mcaMean = mca->GetBinCenter(mca->GetMaximumBin());
  Float_t mcaStdDev = mca->GetRMS();
  fitFunc->SetParameter(pAmpl, mcaPeak);
  fitFunc->SetParameter(pMean, mcaMean);
  fitFunc->SetParameter(pStdDev, mcaStdDev*0.7);
  //Int_t pDCOffset = 3;
  //fitFunc->SetParameter(pDCOffset, 0.0);

  Float_t nsigma = 2.0;
  fitFunc->SetRange(mcaMean-nsigma*mcaStdDev, mcaMean+nsigma*mcaStdDev);
  //Float_t minVal = mca->GetBinLowEdge(1);
  //Float_t maxVal = mca->GetBinLowEdge(mca->GetNbinsX());
  //fitFunc->SetRange(minVal, maxVal);
  //cout << "minVal, maxVal = " << minVal << ", " << maxVal << endl;
  cout << "ampl   = " << fitFunc->GetParameter(pAmpl) << endl;
  cout << "mean   = " << fitFunc->GetParameter(pMean) << endl;
  cout << "stddev = " << fitFunc->GetParameter(pStdDev) << endl;
  //cout << "dc     = " << fitFunc->GetParameter(pDCOffset) << endl;
  //mca->Fit("fitFunc", "L"); // Log likelihood
  mca->Fit("fitFunc", "R");

  return fitFunc;
}

Float_t calcMCAMean(Int_t runnum) {

  TH1F *mca = makeMCA(runnum);
  
  TF1 *f1 = fitToMCA(mca);
  Float_t mean = f1->GetParameter(pMean);
  cout << "mean = " << mean << endl;
  return mean;
}

TString fileOfRun(Int_t runnum, TString fileroot) {
  stringstream ss;
  ss << setfill('0') << setw(6) << runnum;
  TString filename = TString(fileroot);
  filename += ss.str();
  filename += ".root";
  return filename;
}

void processListOfFiles(TString filename) {
  // filename is the name of an ASCII file containing
  // root filenames, one per line, to be processed
  TObjArray *filelist = getFilenames(filename);
  TIter next(filelist);
  while (TObjString* rootFile=(TObjString*)next()) {
    cout << rootFile->GetString() << endl;
    processFile(rootFile->GetString());
  }
}

void processFileRange(TString fileroot, Int_t min, Int_t max) {
  
  for (Int_t ii=min; ii<=max; ii++) {
    cout << "ii = " << ii << "  ";
    TString filename = fileOfRun(ii);
    cout << filename << endl;
    processFile(filename);
  }
}

void processFile(TString infile) {
  
  //if (!fileIsValid(infile)) {
  //  cout << "invalid filename:  " << infile << endl;
  //  return;
  //}

  TString outfile = TString(infile);
  outfile.ReplaceAll(".root", "_red.root");
  cout << "Output: " << outfile << endl;

  // loop over all events in the raw data file
  // compute parameters of interest (ntrigs, pulse heights, etc.)
  DmtpcDataset *d = new DmtpcDataset();
  d->openRootFile(infile);

  Int_t nevents = d->tree()->GetEntries();
  cout << "nevents = " << nevents << endl;

  Float_t pulseHeight;
  TFile *file    = new TFile(outfile, "RECREATE", "reduced cv1 data");
  TTree *outtree = new TTree("outtree", "outtree");
  outtree->Branch("pulseHeight", &pulseHeight, "pulseHeight/F");

  Int_t nWfInEvent, nTrigInEvent;
  Int_t nTrigTotal=0;
  for (Int_t iev=0; iev<nevents; iev++) {
    d->getEvent(iev);
    // nWfInEvent is the number of stored waveforms
    // But, since there are 2 scope channels, this is 2x the number of triggers
    nWfInEvent   = d->event()->scopeData()->GetEntries();
    nTrigInEvent = nWfInEvent/nChan;
    // the first  nWfInEvent/2 entries in scopeData() are from channel 1
    // the second nWfInEvent/2 entries in scopeData() are from channel 2
    // we don't care about channel 2 (nothing was connected)
    for (Int_t itrg=0; itrg<nTrigInEvent; itrg++) {
      nTrigTotal++;
      MaxCamWaveformTools wfA(d->event()->scopeData(itrg));
      pulseHeight = wfA.getTroughDepth();
      outtree->Fill();
    }
  }
  cout << "nTrigTotal = " << nTrigTotal << endl;

  outtree->Write();
  file->Close();
}

bool fileIsValid(TString filename) {
  cout << "Input:  " << filename << endl;
  return true;
}

TObjArray * getFilenames(TString listname, Bool_t verbose) {

  TObjArray *rootList = new TObjArray();
  ifstream ifs(listname.Data());
  string line;
  while (getline(ifs, line)) {
    //if (verbose) cout << "[" << line << "]" << endl;
    TObjString *filename = new TObjString(line.c_str());
    rootList->Add(filename);
  }

  if (verbose) {
    TIter next(rootList);
    while (TObjString* myRootFile=(TObjString*)next()) {
      cout << myRootFile->GetString() << endl;
    }
  }

  return rootList;
}

void plotResolution(TString filename) {
  TObjArray* filelist = getFilenames(filename);
  TIter next(filelist);
  Int_t nfiles = filelist->GetEntries();
  TArrayF means(nfiles);
  TArrayF stdevs(nfiles);
  TArrayF ratios(nfiles);
  Int_t ii = 0;
  while (TObjString* rootFile=(TObjString*)next()) {
    TH1F *mca = makeMCA(rootFile->GetString());
    TF1 *f1   = fitToMCA(mca);
    means.AddAt(f1->GetParameter(1), ii);
    stdevs.AddAt(f1->GetParameter(2), ii);
    ratios.AddAt(f1->GetParameter(2)/f1->GetParameter(1), ii);

    ii++;
  }

  cout << "mean, stdev, ratio = " << endl;
  for (Int_t ii=0; ii<nfiles; ii++) {
    cout << means.At(ii) << "   " << stdevs.At(ii) << "    " << ratios.At(ii) << endl;
  }

}
