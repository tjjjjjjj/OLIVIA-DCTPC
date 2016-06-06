#include <iostream>
#include "TTimeStamp.h"
#include "TDatime.h"
#include "TGraph.h"

using namespace std;

//year, month, day are correct, but there looks to be a timezone issue...

void timetest() {

  TTimeStamp t0;
  TDatime td0;
  const Int_t NN=100;
  Double_t times[NN];
  Double_t yy[NN];
  for (Int_t ii=0; ii<NN; ii++) {
    //TTimeStamp tt(UInt_t(t0.AsDouble()+ii*1000));
    TTimeStamp tt(UInt_t(t0.AsDouble()+ii));
    TDatime td(td0.Convert()+ii);
    //times[ii] = tt.AsDouble();
    times[ii] = Double_t(td.Convert());
    yy[ii]    = Double_t(ii*ii);
  }

  TGraph *gr = new TGraph(NN, times, yy);
  gr->Draw("AP");
  gr->GetXaxis()->SetTimeDisplay(1);
  TTimeStamp epoch1970(1970, 1, 1, 0, 0, 0);
  gr->GetXaxis()->SetTimeOffset(UInt_t(epoch1970.AsDouble()));
  gr->GetXaxis()->SetTimeFormat("#splitline{%d/%m/%y}{%H:%M:%S}");

  return;
}
