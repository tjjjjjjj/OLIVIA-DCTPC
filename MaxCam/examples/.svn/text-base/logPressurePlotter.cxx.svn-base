#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TH1F.h"
using namespace std;

void logPressurePlotter(Int_t refreshRate=1000) {
  // refresh rate is time between plot updates in milliseconds

  TFile f("pressure.root");
  TTree *tree = (TTree*)f.Get("tree");
  TCanvas c1;
  while(1) {
    //tree->Draw("prsCDG", "prsCDG>0");
    //tree->Draw("prsCDG:time>>hnew", "prsCDG>0");
    tree->Draw("prsCDG:time", "prsCDG>0");
    //TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
    //htemp->GetXaxis()->SetNdivisions(10);
    c1.Update();
    gSystem->Sleep(refreshRate); // sleep one second
    tree->Refresh();
  }

}
