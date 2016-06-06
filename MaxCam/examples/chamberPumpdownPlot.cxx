// This is almost there.

// We are not dealing with timezones correctly though.
// The data is plotted with a 5 hour difference, but in the wrong direction.
// e.g. data taken at 11am local time (Boston) is displayed as 06:00 (am)
// even though UTC would be 16:00...
// Not sure if this problem is in data storage or display...

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTimeStamp.h"

using namespace std;

void plot(bool refreshRate_ms=0) {

  // refreshRate_ms is the time interval between updates of the 
  // plotting window.
  // if 0 then only plot once and return.
  
  //Int_t refreshRate=1000; // milliseconds
  //Int_t refreshRate=5000; // milliseconds
  
  TFile f("chamberPumpdown.root");
  TTree *tree = (TTree*)f.Get("tree");
  TCanvas c1;
  TTimeStamp epoch1970(1970, 1, 1, 0, 0, 0);
  TTimeStamp epoch1970A(UInt_t(epoch1970.AsDouble()), kFALSE);
  TTimeStamp epoch1970B(UInt_t(epoch1970.AsDouble()), kTRUE);
  cout << epoch1970.AsString("s") << endl;
  cout << epoch1970A.AsString("s") << endl;
  cout << epoch1970B.AsString("s") << endl;
    
  TTimeStamp *now = new TTimeStamp();
  while(1) {
    now->Set();
    cout << now->AsDouble() << endl;

    TString cut("");
    cut += "time>(";
    cut += now->AsDouble();
    cut += "-7200.)";
    cout << "cut = " << cut << endl;
    tree->Draw("prsBPG:time", cut);
    TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetTimeDisplay(1);
    htemp->GetXaxis()->SetTimeFormat("#splitline{%d/%m/%y}{%H:%M:%S}");
    //htemp->GetXaxis()->SetTimeOffset(0);
    htemp->GetXaxis()->SetTimeOffset(UInt_t(epoch1970.AsDouble()));
    htemp->SetMarkerStyle(6);
    //htemp->GetXaxis()->SetTimeOffset(UInt_t(epoch1970B.AsDouble()));

    c1.Update();
    if (refreshRate_ms > 0) {
      gSystem->Sleep(refreshRate_ms);
      tree->Refresh();
    } else {
      break;
    }
  }
}
