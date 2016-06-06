#include "../DmtpcDataset.hh"
#include "../DmtpcEvent.hh"
#include "../MaxCamSRIM.hh"
#include "../MaxCamConfig.hh"

#include "TChain.h"
#include "TGraph.h"

#include <iostream>
using std::cout;
using std::endl;

#include <vector>

DmtpcDataset *d;

void addRun(int i) {
    TString fname="data/dmtpc_run";
    if (i<100) fname+="000";
    else if (i<1000) fname+="00";
    else if (i<10000) fname+="0";
    fname += i;
    fname += ".root";
    d->chain()->Add( fname );
}

void init() {

  d= new DmtpcDataset;
  d->openRootFile("data/dmtpc_run00801.root");
  for (int i=802; i<=969; i++) {
      addRun(i);
  }

  
  d->getEvent(0);
  TString T0="%b/%d%F";
  T0 += d->event()->timeStamp()->GetYear();
  T0 += "-";
  T0 += d->event()->timeStamp()->GetMonth();
  T0 += "-";
  T0 += d->event()->timeStamp()->GetDay();
  T0 += " ";
  T0 += d->event()->timeStamp()->GetHour();
  T0 +=":";
  T0 += d->event()->timeStamp()->GetMinute();
  T0 +=":";
  T0 += d->event()->timeStamp()->GetSecond();
  
  float pressure = d->event()->experimentConfig("pressure")->currentValue;
  float temperature = d->event()->experimentConfig("temp0")->currentValue+273.15;
  float rho=MaxCamSRIM::density(pressure, 88, temperature); // g/cm3
  float volume = 14.5*14.5*19.5 + 16.5*16.5*19.5; // cm3
  float m=volume*rho; // in g
  cout << "Pressure .... " << pressure << " Torr" << endl;
  cout << "Volume ...... " << volume << " cm3" << endl;
  cout << "Density ..... " << rho << " g/cm3" << endl;
  cout << "Mass ........ " << m << "g"<<endl;

  vector<float> x, y;
  int n=0;


  time_t t0=0;
  time_t event_time=0;
  float tot=0;
  int eventStep=1000;
  float totExpo=0;
  for (int i=0; i<d->chain()->GetEntries(); i+=eventStep) {
    d->getEvent(i);
    event_time=d->event()->timeStamp()->Convert();
    if (i==0) t0=event_time;
    event_time = event_time - t0;
    float eventExposure  = d->event()->ccdConfig(0)->exposureTime*1e-3; // s
    float exposureTime = eventExposure * eventStep /  3600.; // h
    if (n) totExpo+=exposureTime;
    tot += exposureTime / 24  * m; // g*day
    x.push_back(event_time);
    y.push_back(tot);
    cout << " run=" << d->event()->runNumber() 
	 << "  T0=" << event_time/3600. << "h"
	 << "  dT=" <<  exposureTime << "h"
	 << "  exposure="<< tot<< "g-day";
    if (n) cout << "   eff="<< totExpo/event_time*3600;
    cout << endl;
    n++;
  }
  TGraph *lumi=new TGraph(x.size()-1,&x[1],&y[1]);
  lumi->SetLineWidth(3);
  lumi->SetLineColor(4);
  lumi->Draw();
  lumi->GetXaxis()->SetTimeDisplay(1);
  lumi->GetXaxis()->SetTimeFormat((const char*)T0);
  lumi->Draw("AL");
  lumi->GetHistogram()->SetYTitle("Exposure (g-day)");

}
