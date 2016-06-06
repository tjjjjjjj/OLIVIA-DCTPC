#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "../MaxCamSRIM.hh"

using namespace std;

void makeplot() {

  MaxCamSRIM srimHe("SRIM_He_in_CF4_100Torr");

  double pMin = 1.0; // torr
  double pMax = 800.0; // torr
  int np = 100;  // number of samples in pressure
  double dp = (pMax-pMin)/(np-1);
  double pp;
  double range;
  TGraph alpha5MeV(np);
  TGraph alpha4MeV(np);
  TGraph alpha6MeV(np);

  for (int ii=0; ii<np; ii++) {
    pp = pMin +(ii*dp);
    //cout << pp << " " << endl;
    srimHe.setPressure(pp);
    range = srimHe.getRangeVsEnergy()->Eval(5000);
    alpha5MeV.SetPoint(ii, pp, range*0.1);

    range = srimHe.getRangeVsEnergy()->Eval(4000);
    alpha4MeV.SetPoint(ii, pp, range*0.1);

    range = srimHe.getRangeVsEnergy()->Eval(6000);
    alpha6MeV.SetPoint(ii, pp, range*0.1);


    //cout << "pp, range (mm) = " << pp << ", " << range << endl;
  }
  TCanvas *c1 = new TCanvas("c1", "example", 200, 10, 700, 500);
  alpha5MeV.Draw("AC");
  alpha5MeV.GetHistogram()->SetXTitle("Pressure [torr]");
  alpha5MeV.GetHistogram()->SetYTitle("Range [cm]");

  alpha6MeV.SetLineStyle(2);
  alpha6MeV.Draw("CSAME");

  alpha4MeV.SetLineStyle(3);
  alpha4MeV.Draw("CSAME");
  
  TLegend *legend = new TLegend(0.60, 0.68, 0.75, 0.8);
  legend->AddEntry(&alpha5MeV, "5 MeV in CF4", "L");
  legend->AddEntry(&alpha4MeV, "4 MeV", "L");
  legend->AddEntry(&alpha6MeV, "6 MeV", "L");
  legend->Draw();

  c1->SetLogy(1);
  c1->SetTitle("Range vs. Pressure for Alphas in CF4");
  c1->Print("test.pdf");

  //Double_t xx, yy;
  //for (Int_t ii=0; ii<10; ii++) {
  //  alpha5MeV.GetPoint(ii, xx, yy);
  //  cout << "xx, yy = " << xx << ", " << yy << endl;
  //}

}


