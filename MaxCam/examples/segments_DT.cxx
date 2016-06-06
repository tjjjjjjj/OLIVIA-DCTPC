//////////////////////////////////////////////////////////////////////////////////
//
//   Recoinstruction of recoil segments for neutron scattering (DT) experiment
//   Data runs 139+146+149
//   Publications: NIM/EPS conf
//
//   - in ROOT:
//   .L segments.cxx++
//   init()
//   segments()
//   - output file is segout.root
//   - make plots using segmentsPlot_DT.cxx
//
//////////////////////////////////////////////////////////////////////////////////


#include "../MaxCamRead.hh"
#include "../MaxCamMC.hh"
#include "../MaxCamTrack.hh"
#include "TCanvas.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TF1.h"
#include "TTree.h"
#include <iostream>
using std::cout;
using std::endl;

bool makeDark=false; // dark plots for presentations

int i=0;
MaxCamRead *ana=0;
TCanvas *c1, *c2;
void init() {
        gSystem->Setenv("MCTABLES","/Users/ddujmic/work/projects/DarkMatter/MaxCam/tables");

        // neutron runs (0, 180deg): 139, 146, 149
        // alpha calibration runs: 154-162, 171
        // neutrons perp to wires (90deg): 185-188
        // neutrons with wall-gap closed (147)
        ana = new MaxCamRead("data/ccdrun_00294.root");
        
        ana->readBiasFrame();
        ana->findHotPixels("hotpixels.dat");
        
        //ana = new MaxCamRead("MonteCarloSamples/mcrun_00054.root");
        //ana = new MaxCamRead("mcrun.root");
        
        c1=new TCanvas("c1","",0,0,700,600);
        c2=new TCanvas("c2","",700,0,750, 600);
        if (makeDark) {
            gStyle->SetLineWidth(2);
            c1->SetFillColor(1);
            c1->SetFillStyle(0);
            c2->SetFillColor(1);
            c2->SetFillStyle(0);
        }
        c2->Divide(3,2);
        c1->cd();

        ana->findWires(2, 0.4, 100, "y"); // 139,146,149 recoil
        ana->srim->fillSrimTable("SRIM_F_in_CF4_100Torr");
        ana->srim->setStopping(200);
        ana->srim->setProjectile(-14.1e3, 0, 0, 1e6); // 14.1MeV neutron
        ana->srim->setRecoil(100, 0, 0, 19e6); // 100keV fluorine

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
        if (ana->getVerbose()) cout << "EVENT="<< i<<endl;
        ana->getRunNumber();
        ana->setEventInfo(i-1);
        if (opt.Contains("p")) { c1->cd(); image(); }

        return true;
        //if (ana->ccdImage()->Integral()>1e6) return false;     
        //cout << "CCD integral="<<   ana->ccdImage()->Integral() << endl;
        //if (ana->ccdImage()->Integral()<1) return false;

        
        TH1F *hY = ana->makeYieldHisto();
        MaxCamTrack* trfit = new MaxCamTrack( ana->ccdImage(), true );
        float threshold = hY->GetMean() + hY->GetRMS()*3;
        trfit->setThreshold( threshold );
        trfit->setWireBinList( ana->wireBinList, "x" ); 
        //trfit->setWireBinList( ana->wireBinList, "y" ); // recoils
        //trfit->checkThresholdAndNeighbors(); trfit->getFitImage()->Draw("colz"); return true;
        trfit->makeTracks(); 
        cout << ": Found # of tracks="<< trfit->nTracks() << endl;
        if (trfit->nTracks()!=1) return false; // recoils
        //if (trfit->nTracks()<1) return false;
        float slope=trfit->getTrack(0)->GetParameter(1);
        if (fabs(slope)>0.1) return false;
        cout << "slope="<<trfit->getTrack(0)->GetParameter(1)<< endl;
        delete hY;
        delete trfit;
        
        return true;
}


TH1* strip;
void dolength(int iwire, int dwire=2, float threshold=0.0, int gapmax=0, TString opt="") {

    //strip = ana->findSegmentLength(iwire, dwire,threshold,gapmax,"xrc");
    strip = ana->findSegmentLength(iwire, dwire, threshold, gapmax,"yrc"); // recoils
    
    if (opt.Contains("p")) {
        TH1F *stripclone=(TH1F*)strip->Clone("stripclone");
        stripclone->SetFillColor(1);
        stripclone->SetFillStyle(1001);
        strip->GetXaxis()->SetRange( 0,0);
        strip->Draw();
        stripclone->Draw("same");
    }
    else  delete strip;
    
}



int analyze(TString opt="") {
    
    c2->cd(1);
    TH1F *hY=ana->makeYieldHisto();
    float threshold = hY->GetMean() + hY->GetRMS()*3;
    delete hY;

    ana->setFitNumber(0);
    for (unsigned int iw=0; iw<ana->wireList.size(); iw++) {
        c2->cd(iw+1); dolength(iw, 2, threshold,  0, opt);
        c2->cd(iw+4); ana->wireMinROI=(int)ana->getWidthL(iw); ana->wireMaxROI=(int)ana->getWidthR(iw);
        if (ana->getWidthL(iw)<ana->getWidthR(iw)) { ana->wireMinROI-=5;  ana->wireMaxROI+=5; }
        ana->wireYield(iw,10,1,opt); 
    }
    fill();
    
    return 0;
}








void segments() {

    
    int ev139[]={10,46,60,62,90,112,121,124,127,131,165,182,185,225,243,256,260,272,305,315,317,319,342,351,356,379,391,
                 404,406,409,426,439,452}; // run 139

    int ev149[]={4,8,9,16,19,36,45,80,86,95,107,206,225,244,252,288,307,309,311,315,322,332,355,410,413,418,419,425,
              427,444,451,487,490,525,538,543,556,572,574,581,586,609,617,623,640,641,642,653,666,675,685,696,711,721,791}; //run 149

    int ev146[]={8,9,13,18,20,21,27,53,69,70,94,146,158,199,231,282,296,301,320,329,330,351,369,374,387,391,401,406,456,
              466,469,479,503,515,522,527,530,555,561,570,611,617,621,638,685,687,689,706,718,730,744,752,757,758,782,796,816,
              843,902,907,940,972,988, 1034,1037,1038,1044,1075,1077,1088,1090,1102,1110,1121,1137,1140,1152,1156,1187,
              1312,1316,1373,1379,1417,1422,1443,1485,1490,1514,1624,1631,1644,1653,1658,1668,1706,1744,1754,1780,1809,
              1823,1857,1904,1906,1917,1927,1942,1944,1949,1951,1987,2051,2063,2098,2123,2168,2174,2181,2208,
              2214,2220,2295,2311};  // run146
    

    
    int *ev=0, nev=0;
    switch (ana->getRunNumber()) {
        case 139: ev=ev139; nev=sizeof(ev139)/sizeof(int); break;
        case 146: ev=ev146; nev=sizeof(ev146)/sizeof(int); break;
        case 149: ev=ev149; nev=sizeof(ev149)/sizeof(int); break;
    }

    
    nev=ana->tree()->GetEntries();
    cout << "Selected events="<<nev<<endl; 
    //cout << "NOT STARTING FORM 1ST EVENT!!!!!!!!" << endl;
    for (int i=0; i<nev; i++) {
        //if (!event( ev[i]) ) continue;
        if (!event(-1 )) continue;
        //event( ev[i] );
        //event(-1);
        //analyze("y");
        analyze("x"); // recoils
    }
    TFile f("segout.root","recreate");
    ana->nt->Write();
    f.Write();
}



/////////////////////////////////////
//
//   NIM plots
//
/////////////////////////////////////
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "../MaxCamImageTools.hh"
void makeFigure6_NIM(int j=0) {

    
    //int col[]={12,13,14,15,16,17,18,19,10};
    //gStyle->SetPalette(sizeof(col)/sizeof(int),col);
    gStyle->SetPadTopMargin(0.05); 
    gStyle->SetPadRightMargin(0.05); 
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetPadLeftMargin(0.1);
    if (!ana) init();

    //c1->SetTheta(28.0919);
    //c1->SetPhi(18.7663);

    
    int run = ana->getRunNumber();
    
    int event139[]={182,315,342};
    int xroi139[]={ 7,51, 39,89, 10,52};
    int yroi139[]={46,54, 6,14,  46,54};

    int event149[]={311,543,572};
    int xroi149[]={ 0,69,   39,89,  13,95};
    int yroi149[]={ 53,61,  53,61,  53,61};

    int *evt=0, *xroi=0, *yroi=0;
    switch (run) {
        case 139: evt=event139; xroi=xroi139; yroi=yroi139; break;
        case 149: evt=event149; xroi=xroi149; yroi=yroi149; break;
        default: assert(0);
    }

    
    TString fname="run00";
    fname+=run;
    fname+="_ev";
    fname+=evt[j]+1;
    fname+="_recoilCandidate.pdf";
    cout << fname << endl;
    event(evt[j],"p");
    ana->ccdImage()->GetXaxis()->SetRange(xroi[2*j],xroi[2*j+1]);
    ana->ccdImage()->GetYaxis()->SetRange(yroi[2*j],yroi[2*j+1]);
    ana->ccdImage()->GetXaxis()->SetNdivisions(505);
    ana->ccdImage()->GetYaxis()->SetNdivisions(505);
    ana->ccdImage()->GetZaxis()->SetNdivisions(505);
    ana->ccdImage()->Draw("lego2 fb");
    if (makeDark) { // e.g. presentations
        ana->ccdImage()->SetLineColor(10);
        ana->ccdImage()->SetLineWidth(3);
        ana->ccdImage()->GetXaxis()->SetLabelColor(10);
        ana->ccdImage()->GetXaxis()->SetTitleColor(10);
        ana->ccdImage()->GetXaxis()->SetAxisColor(10);
        ana->ccdImage()->GetYaxis()->SetLabelColor(10);
        ana->ccdImage()->GetYaxis()->SetTitleColor(10);
        ana->ccdImage()->GetYaxis()->SetAxisColor(10);
        ana->ccdImage()->GetZaxis()->SetLabelColor(10);
        ana->ccdImage()->GetZaxis()->SetTitleColor(10);
        ana->ccdImage()->GetZaxis()->SetAxisColor(10);
        c1->Print("~/Desktop/c1.gif");
    }
    c1->Print(fname);
}

void makeFigure2_NIM() {

    //int col[]={12,13,14,15,16,17,18,19,10};
    int col[]={10, 19,18,17,16,15,14,13,12, 1};
    gStyle->SetPalette(10,col);
    

    ana = new MaxCamRead("data/ccdrun_00171.root");
    
    ana->readBiasFrame();
    ana->findHotPixels("hotpixels.dat");
    
    TH2F* h=ana->accumulatePressure(280);
    h->SetBinContent(34,19,0);
    TH2F* hmm=MaxCamImageTools::resizeImage(h, 0, 768/22.,  0, 512/22.);
    hmm->Draw("colz");
    hmm->SetXTitle("x (mm)");
    hmm->SetYTitle("y (mm)");


    TLine *l=new TLine(6,20, 2.6,22.5);
    l->SetLineWidth(10);
    l->Draw();
    TLatex *lat=new TLatex(4.,21.4,"^{241}Am");
    lat->SetTextSize(0.05);
    lat->Draw();

    c1->Print("run00171_allEvents_mm_bw.pdf");
}
