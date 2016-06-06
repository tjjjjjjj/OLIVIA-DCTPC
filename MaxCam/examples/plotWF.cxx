#include "TNtuple.h"
#include "TString.h"
#include "TCanvas.h"
#include "../DmtpcDataset.hh"
#include "../MaxCamWaveformTools.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TStyle.h"

#include <iostream>
using std::cout;
using std::endl;

TNtuple *nt=0;
DmtpcDataset *d=0;
TCanvas *c=0;

void init(TString fname="test.root", TString opt="") {

  if (!d) {
    d= new DmtpcDataset;
    d->openRootFile(fname);
  }

  opt.ToLower();

  TFile resfile("veto.root","RECREATE");
  nt=new TNtuple("nt", "", "run:event:maxi0:maxi1:time0:time1:rise0:rise1:ccd:drift:anode:nwf:mini0:mini1:diff0:diff1");


  c= new TCanvas;

  TStopwatch sw;
  sw.Start();

  int nev=d->tree()->GetEntries();
  for (int i=0; i<nev; i++) {
      if (i%10==0) cout << "event " << i << endl;
    d->tree()->GetEvent(i);
    int nwf=d->event()->scopeData()->GetEntries()/2;

    float drift=d->event()->experimentConfig("driftHV") ? d->event()->experimentConfig("driftHV")->setValue : 0;
    float anode=d->event()->experimentConfig("anodeHV") ? d->event()->experimentConfig("anodeHV")->setValue : 0;
    
    for (int j=0; j<nwf; j++) {

        TH1F *chA=d->event()->scopeData(j);
        TH1F *chB=d->event()->scopeData(j+nwf);
        
        MaxCamWaveformTools wfA(chA); //wfA.print();
        MaxCamWaveformTools wfB(chB); //wfB.print();

        float var[]={d->event()->runNumber(),
                     d->event()->eventNumber(),
                     wfA.getPeakHeight(),
                     wfB.getPeakHeight(),
                     wfA.getPeakHeightTime(),
                     wfB.getPeakHeightTime(),
                     wfA.getRiseTime(),
                     wfB.getRiseTime(),
                     d->event()->ccdData(0)->Integral(),
                     drift, anode,
                     d->event()->scopeData()->GetEntries(),
                     wfA.getWfMinimum(),
                     wfB.getWfMinimum(),
                     wfA.getPeakAveragePosTime()-wfA.getPeakAverageNegTime(),
                     wfB.getPeakAveragePosTime()-wfB.getPeakAverageNegTime(),                  
        };
        
        nt->Fill( var ); 
        
      

    }
    
  }
  sw.Print();

  resfile.Write();
  resfile.Close();
  delete d; d=0;
}

#include "TChain.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"

void peak() {

    TH1F h1("h1","",128,0,2);
    TH1F h2("h2","",128,0,2);
    TH1F hdiff("hdiff","",128,0,2);

    TChain nt("nt");
    nt.Add("veto_run01463.root");
    nt.Add("veto_run01464.root");

    nt.Draw("maxi0>>h1","rise0<5e-6&&rise0>0&&run==1463");
    nt.Draw("maxi0>>h2","rise0<5e-6&&rise0>0&&run==1464");
    hdiff.Add(&h1, &h2, 1, -1);
    hdiff.DrawCopy();
}



void rate() {

    TChain nt("nt");
    nt.Add("veto_run01463.root");
    TString cut0="rise0<5e-6&&rise0>0&&maxi0<2.4";

    vector<float> x,y;
    for (float ecut=0; ecut<2; ecut+=0.1) {
        TString cut=cut0;
        cut+="&&maxi1>";
        cut+=ecut;
        x.push_back( ecut/1.3*350 ); // calibration from am peak in veto
        y.push_back(  nt.Draw("1", cut)*1e-3 ); // 1000 events total
    }

    gPad->Clear();
    TGraph *g=new TGraph(x.size(), &x[0], &y[0]);
    g->SetMaximum(1);
    g->Draw("AP");
}


void calcRiseTime(TTree *nt, int irun, float &mean, float &rms) {

    TString cut="maxi0>0.05&&maxi0<0.99&&run==";
    cut += irun;

    TH1F ht("ht","",100,0.0, 0.5e-6);
    nt->Draw("rise0>>ht", cut);

    
    mean=ht.GetMean();
    rms=ht.GetRMS();
    ht.Fit("gaus","LL");
    mean=ht.GetFunction("gaus")->GetParameter(1);
    rms=ht.GetFunction("gaus")->GetParameter(2);
}


void plotRiseVsAngle() {

    TChain nt("nt");
    //nt.Add("veto_run01468.root");
    //nt.Add("veto_run01469.root");
    //nt.Add("veto_run01470.root");
    //nt.Add("veto_run01474.root");  
    //float irun[]={1468, 1469, 1470,  1474};
    //float angle[]={30, 45, 60, 0};

    nt.Add("veto_run01479.root");
    nt.Add("veto_run01480.root");
    float irun[]={1479, 1480};
    float deltaz[]={0, 7};
    vector<float> y, yrms, x, xerr;

    int n=sizeof(irun)/sizeof(int);
    float mean, rms;
    for (int i=0; i<n; i++) {
        x.push_back( deltaz[i] );
        xerr.push_back( 1 );

        calcRiseTime( &nt, int(irun[i]), mean, rms);
        y.push_back( mean*1e9/0.8 );
        yrms.push_back( rms*1e9/0.8 );
        
    }

    TGraphErrors *g=new TGraphErrors(x.size(), &x[0],&y[0], &xerr[0], &yrms[0]);
    gPad->Clear();
    g->Draw("AP");
}


void PulseVsRise() {
    gStyle->SetOptTitle();
    TChain nt("nt");
    nt.Add("veto_run01532.root");  // DT
    nt.Add("veto_run01529.root");  // Cs137
    nt.Add("veto_run01530.root");  // none
    
    TCanvas *c=new TCanvas("c");
    c->Divide(2,2);
    
    TH2F *frame=new TH2F("frame","", 100,0, 1.5, 128,0,1);
    frame->SetXTitle("Rise time (#mus)");
    frame->SetYTitle("Pulse height (V)");


    c->cd(1);
    nt.Draw("maxi0:rise0*1e6>>frame", "maxi0>0.05&&run==1530");
    frame->Scale(1./frame->Integral());
    frame->SetTitle("run 1530 (bkg)");
    frame->DrawCopy("colz");
    
    c->cd(2);
    nt.Draw("maxi0:rise0*1e6>>frame", "maxi0>0.05&&run==1532");
    frame->Scale(1./frame->Integral());
    frame->SetTitle("run 1532 (DT)");
    frame->DrawCopy("colz");

    c->cd(3);
    nt.Draw("maxi0:rise0*1e6>>frame", "maxi0>0.05&&run==1529");
    frame->Scale(1./frame->Integral());
    frame->SetTitle("run 1529 (Cs137)");
    frame->DrawCopy("colz");


    
}
