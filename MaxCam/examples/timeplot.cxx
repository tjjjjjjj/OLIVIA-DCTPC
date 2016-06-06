#include <iostream>
#include "TDatime.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TGraph.h"
#include "TCanvas.h"
using namespace std;

void timeplot() {

  TTimeStamp time;
  //TDatime time;

  const int NN = 100;
  Double_t timedbl[NN];
  Double_t yy[NN];

  TTimeStamp epoch1995(1995, 1, 1, 0, 0, 0);
  TTimeStamp epoch1970(1970, 1, 1, 0, 0, 0);
  Double_t d95_70 = epoch1995.AsDouble()-epoch1970.AsDouble();

  cout << epoch1995.AsString() << endl;
  cout << epoch1970.AsString() << endl;

  //TDatime t0(time);
  TTimeStamp t0;
  for (int ii=0; ii<NN; ii++) {
    time.Set();
    //timedbl[ii] = time.Convert();
    timedbl[ii] = time.AsDouble() - d95_70;  // return time since 1995
    yy[ii] = 1.0*ii*ii;
    gSystem->Sleep(50);
  }

  cout << time.AsString() << endl;

  TCanvas *c1 = new TCanvas();
  TGraph  *gr = new TGraph(NN, timedbl, yy);
  gr->Draw("AP");
  gr->GetXaxis()->SetTimeDisplay(1);
  
  TString tfmt = "#splitline{%d/%m/%y}{%H:%M:%S %p}";
  //cout << "t0.AsSQLString() = " << t0.AsSQLString() << endl;
  //cout << "t0.Convert() = " << t0.Convert() << endl;
  cout << "t0.AsString() = " << t0.AsString() << endl;
  cout << "UInt_t(t0.AsDouble()) = " << UInt_t(t0.AsDouble()) << endl;
  gr->GetXaxis()->SetTimeFormat(tfmt.Data());

  UInt_t offset = UInt_t(t0.AsDouble() - epoch1995.AsDouble());
  cout << "offset = " << offset << endl;
  //gr->GetXaxis()->SetTimeOffset(offset);
  //gr->GetXaxis()->SetTimeOffset(UInt_t(epoch1995.AsDouble()));
  
  c1->Update();

  return;
}
