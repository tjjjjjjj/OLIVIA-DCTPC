#include "../MaxCamRead.hh"
#include "../MaxCamMC.hh"
#include "../MaxCamTrack.hh"
#include "TCanvas.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TF1.h"
#include "TTree.h"
#include <iostream>
using std::cout;
using std::endl;


int i=0;
MaxCamRead *ana;
TCanvas *c1, *c2;
void init() {
        gSystem->Setenv("MCTABLES","../tables");

        ana = new MaxCamRead("data/ccdrun_00231.root");
        ana->readBiasFrame();
        ana->findHotPixels("hotpixels.dat");
        //ana->setDay0(9);
        
        c1=new TCanvas("c1","",0,0,700,600);
        c2=new TCanvas("c2","",700,0,750, 600);
        c2->Divide(2,2); 
        c1->cd();

        ana->findWires(4, 0.3, 400, "x");
}

void fill() { ana->nt->Fill(); }

void image() { 
        TH2F *hall=ana->ccdImage(); 
        hall->SetMaximum(500); 
        hall->Draw("colz"); 
}


bool event(int ii=-1, TString opt="") {
        if (ii>0) i=ii;
        ana->getEvent(i++); 
        cout << "EVENT="<< i<<endl;
        ana->getRunNumber();
        ana->setEventInfo(i-1);
        if (opt.Contains("p")) { c1->cd(); image(); }

        //cout << "CCD integral="<<   ana->ccdImage()->Integral() << endl;
        if (ana->ccdImage()->Integral()>1e6) return false;     
        if (ana->ccdImage()->Integral()<1) return false;     

        TH1F *hY = ana->makeYieldHisto();
        MaxCamTrack* trfit = new MaxCamTrack( ana->ccdImage() );
        float threshold = hY->GetMean() + hY->GetRMS()*3;
        cout << "THRRESHOLD="<<threshold<<endl;
        trfit->setThreshold( threshold );
        trfit->setWireBinList( ana->wireBinList, "x" ); 
        trfit->makeTracks(); 
        //cout << ": Found # of tracks="<< trfit->nTracks() << endl;
        if (trfit->nTracks()!=1) return false; // recoils
	//if (trfit->nTracks()<1) return false;
        float slope=trfit->getTrack(0)->GetParameter(1);
        if (fabs(slope)>0.1) return false;
        cout << "slope="<<trfit->getTrack(0)->GetParameter(1)<< endl;
        delete hY;

        return true;
}


TH1* strip;
void dolength(int iwire, int dwire=2, float threshold=0.0, int gapmax=0, TString opt="") {

    strip = ana->findSegmentLength(iwire, dwire,threshold,gapmax,"xrc");
    
    if (opt.Contains("p")) {
        TH1F *stripclone=(TH1F*)strip->Clone("stripclone");
        stripclone->SetFillColor(1);
        stripclone->SetFillStyle(1001);
        strip->GetXaxis()->SetRange( 0,0);
        strip->Draw();
        stripclone->Draw("same");
    }
    else delete strip;
}



int analyze(TString opt="") {
    
    c2->cd(1);
    TH1F *hY=ana->makeYieldHisto();
    float threshold = hY->GetMean() + hY->GetRMS()*3;
    delete hY;

    for (int i=1; i<2; i++) {
        ana->setFitNumber(i);
        ana->getStoppingFunction()->SetParameters(170, 730, 1.7, -170, 10, pow(-1,i), 5./110);
        c2->cd(1); dolength(0, 2, threshold,  0, opt); 
        c2->cd(2); dolength(1, 2, threshold, 0, opt);  //if (!i) ana->copyShapeToAntishape(1);
        c2->cd(3); ana->wireMinROI=ana->getWidthL(0)-15; ana->wireMaxROI=ana->getWidthR(0)+15;  ana->wireYield(0,10,1,opt); 
        c2->cd(4); ana->wireMinROI=ana->getWidthL(1)-15; ana->wireMaxROI=ana->getWidthR(1)+15;  ana->wireYield(1,10,1,opt);
        dolength(2, 2, threshold, 0, opt);
        ana->wireMinROI=ana->getWidthL(2)-15; ana->wireMaxROI=ana->getWidthR(2)+15;  ana->wireYield(2,10,1,opt);
        //c2->cd(1);ana->wireYield(0,3,1,opt);
        //c2->cd(2);ana->wireYield(1,3,1,opt);
        //c2->cd(3);ana->wireYield(2,3,1,opt);
    }       
    fill();
    
    return 0;
}








void segments() {

    int nev=ana->tree()->GetEntries();

    for (int i=0; i<nev; i++) {
        if (!event(-1 )) continue;
        analyze("y");
        //analyze("x"); // recoils
    }
    TFile f("segout.root","recreate");
    ana->nt->Write();
    f.Write();
}


