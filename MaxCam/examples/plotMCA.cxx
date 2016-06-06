
#include "TH1F.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "TFile.h"
#include "TF1.h"
#include "TChain.h"


#include "../DmtpcDataset.hh"
#include "../DmtpcEvent.hh"
#include "../MaxCamChannel.hh"

#include "time.h"
#include <iostream>
using std::cout;
using std::endl;

TH1F *hwall=0;
DmtpcDataset *d=0;
TNtuple *ntWF=0, *ntEvent=0;
TFile *fout;

int groupEvents=10;



void plot(TString fname="test.root", TString what="mca") {

    time_t event_time=0, t0=0;


    if (!d) {
        d= new DmtpcDataset;
        d->openRootFile(fname);
        cout << "events="<<d->tree()->GetEntries()<<endl;
    }
    TFile* fout = new TFile("plotMCA.root","RECREATE");
    ntWF=(TNtuple*)fout->Get("ntWF");
    if (!ntWF) {
        ntWF=new TNtuple("ntWF","","peak:min:rms:skew:rise:fall:time:run");
        ntEvent= new TNtuple("ntEvent","", "multi:peak:time:overflow:press:anode:drift");
    }

  what.ToLower();

  TH1F hevent("hevent","", 128,    0-0.039,   10-0.039);
 
  int n=d->tree()->GetEntries();
  for (int i=0; i<n; i++) {

      // get event
      d->tree()->GetEvent(i);
      int nwf=d->event()->scopeData()->GetEntries();
      int nth=0;
      event_time=d->event()->timeStamp()->Convert();
      if (i==0) t0=event_time;
      event_time = event_time - t0;

      // loop over waveforms
      for (int j=0; j<nwf; j++) {
          if (what=="threshold"&& d->event()->scopeData(j)->GetMaximum()>2.4 &&  d->event()->scopeData(j)->GetMaximum()<9.9) nth++;

          TH1F *waveform=d->event()->scopeData(j);
          
          // get rise/fall time
          waveform->GetXaxis()->SetRange(800,1300);
          //waveform->Fit("gaus", "Q");
          int kmax=waveform->GetMaximumBin();
          int tlow=0, thigh=0;
          for (int k=kmax; k>800&&0; k--) {
              if ( waveform->GetBinContent(k)/waveform->GetMaximum()<0.9 && thigh==0) thigh=k;
              if ( waveform->GetBinContent(k)/waveform->GetMaximum()<0.1 && tlow==0) { tlow=k; break; }
          }
          float rise=(waveform->GetBinCenter(thigh) - waveform->GetBinCenter(tlow))*1e6;
          tlow=thigh=0;
          for (int k=kmax; k<1300&&0; k++) {
              if ( waveform->GetBinContent(k)/waveform->GetMaximum()<0.9 && tlow==0) tlow=k;
              if ( waveform->GetBinContent(k)/waveform->GetMaximum()<0.1 && thigh==0) { thigh=k; break; }
          }
          float fall=(waveform->GetBinCenter(thigh) - waveform->GetBinCenter(tlow))*1e6;

          if (waveform->GetMaximum()-waveform->GetMinimum()>1.25)
              hevent.Fill(waveform->GetMaximum());
          
          // fill waveform ntuple
              ntWF->Fill( waveform->GetMaximum(), 
                          waveform->GetMinimum(),
                          waveform->GetRMS()*1e6, 
                          waveform->GetSkewness()*1e6, 
                          rise,
                          fall,
                          waveform->GetBinContent(0),
                          d->event()->runNumber()
                          );
      }
      
      // fill event ntuple

      if ((i+1)%groupEvents==0) {

          int nmax=hevent.GetMaximumBin();
          float avgpeak=0, peaknorm=0;
          for (int imax=nmax-5; imax<nmax+5; imax++) {
              avgpeak+= hevent.GetBinCenter(imax)*hevent.GetBinContent(imax);
              peaknorm+= hevent.GetBinContent(imax);
          }
          avgpeak/=peaknorm;
          hevent.GetXaxis()->SetRange(nmax-5, nmax+5);
          hevent.Fit("gaus", "Q");
          cout << avgpeak << "  " << hevent.GetFunction("gaus")->GetParameter(1) << endl;
          //hevent.Draw(); c1->Update(); getchar();
          avgpeak = hevent.GetFunction("gaus")->GetParameter(1);

          float press=d->event()->experimentConfig("pressure")->currentValue;
          float anode=d->event()->experimentConfig("anodeHV")->currentValue;
          float drift=d->event()->experimentConfig("driftHV")->currentValue;
          
          ntEvent->Fill( d->event()->scopeData()->GetEntries(),
                         avgpeak,
                         event_time/3600.,
                         hevent.GetBinContent( hevent.GetNbinsX()),
                         press,
                         anode,
                         drift
                         );
          hevent.Reset();
      }
  }
  

  // output results to file
  
  ntWF->Write();
  ntEvent->Write();
  fout->Write();
  //fout->Close();
  //delete fout;
  
  delete d; d=0;
}


void all() {

    d= new DmtpcDataset;
    d->openRootFile("data/dmtpc_run00792.root");
    d->chain()->Add("data/dmtpc_run00793.root");
    d->chain()->Add("data/dmtpc_run00794.root");
    d->chain()->Add("data/dmtpc_run00795.root");
    d->chain()->Add("data/dmtpc_run00796.root");
    plot();
}




/*

MaxCamSRIM *fluor=0;
TF1 *fq;
//double calibFactor=5.9/2.1; // keV/V
double calibFactor=26/2.37; // keV/V
double fquench(double *var, double *par) {
    
    double Eee=var[0];
    
    double E=Eee, Eold=-1e10;
    if (!fluor) fluor=new MaxCamSRIM("SRIM_F_in_CF4_100Torr");

    return Eee*fluor->getStoppingVsEnergy(0)->Eval(Eee)/fluor->getStoppingVsEnergy()->Eval(Eee);
   
    while (fabs(E-Eold)>0.01) {
        Eold=E;
        E = Eee*fluor->getStoppingVsEnergy()->Eval(E)/fluor->getStoppingVsEnergy(0)->Eval(E);
    }
    return E;
}

void calibration() {
    TGaxis *kevee = new TGaxis(0,2e4, 10,2e4,  0, 10*calibFactor, 505, "-");
    kevee->SetName("kevee");
    kevee->Draw();
    gPad->SetLogy();
    
    fq = new TF1("fq",fquench, 0, 150, 1);
    TGaxis *kevf = new TGaxis(0,8e4, 10,8e4,  "fq", 505, "-");
    kevf->SetName("kevf");
    kevf->Draw();
}
*/

