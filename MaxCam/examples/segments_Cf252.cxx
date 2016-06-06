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
        gSystem->Setenv("MCTABLES","/Users/ddujmic/work/projects/DarkMatter/MaxCam/tables");

        ana = new MaxCamRead("data/ccdrun_00251.root");
        ana->addRun("data/ccdrun_00252.root");
        
        ana->readBiasFrame();
        ana->findHotPixels("hotpixels.dat");
        ana->setDay0(26);
        
        //ana = new MaxCamRead("MonteCarloSamples/mcrun_00048.root");
        //ana->addRun("MonteCarloSamples/mcrun_00051.root");
        //ana = new MaxCamRead("mcrun.root");
        
        c1=new TCanvas("c1","",0,0,700,600);
        c2=new TCanvas("c2","",700,0,750, 600);
        c2->Divide(3,2); 
        c1->cd();
        //ana->getProfileFunction()->SetParLimits(3, 6.8, 6.8);
        //ana->getProfileFunction()->SetParameter(3, 6.8);
        //ana->findWires(4, 0.4, 600, "x"); 
        ana->findWires(2, 0.4, 100, "y"); // 139,146,149 recoil
        ana->srim->fillSrimTable("SRIM_F_in_CF4_100Torr");
        ana->srim->setStopping(250);
        ana->getStoppingFunction()->SetParLimits( 0, 0,   640);
        ana->getStoppingFunction()->SetParLimits( 1, 650, 650);
        ana->getStoppingFunction()->SetParLimits( 2, 0,   100);
        ana->getStoppingFunction()->SetParLimits( 4, 10,   10);
        ana->getStoppingFunction()->SetParLimits( 5, -1,   -1);
        ana->getStoppingFunction()->SetParLimits( 6, 0.0454,  0.0454);

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

        //return true;
        //if (ana->ccdImage()->Integral()>1e6) return false;     
        //cout << "CCD integral="<<   ana->ccdImage()->Integral() << endl;
        //if (ana->ccdImage()->Integral()<1) return false;

        
        TH1F *hY = ana->makeYieldHisto();
        MaxCamTrack* trfit = new MaxCamTrack( ana->ccdImage(), true );
        float threshold = hY->GetMean() + hY->GetRMS()*3;
        //cout << "THRRESHOLD="<<threshold<<endl;
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
    else delete strip;
}



int analyze(TString opt="") {
    
    c2->cd(1);
    TH1F *hY=ana->makeYieldHisto();
    float threshold = hY->GetMean() + hY->GetRMS()*3;
    delete hY;

    ana->setFitNumber(0);
    for (unsigned int iw=0; iw<ana->wireList.size(); iw++) {
        c2->cd(iw+1); dolength(iw, 2, threshold,  0, opt);
        c2->cd(iw+4); ana->wireMinROI=(int)ana->getWidthL(iw)-5; ana->wireMaxROI=(int)ana->getWidthR(iw)+5;  ana->wireYield(iw,10,1,opt); 
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
    
    int ev247[]={73, 149,237,298,338,563,566,870,922,1070,1230,1379,1380};

    
    int *ev=0, nev=0;
    switch (ana->getRunNumber()) {
        case 139: ev=ev139; nev=sizeof(ev139)/sizeof(int); break;
        case 146: ev=ev146; nev=sizeof(ev146)/sizeof(int); break;
        case 149: ev=ev149; nev=sizeof(ev149)/sizeof(int); break;
        case 247: ev=ev247; nev=sizeof(ev247)/sizeof(int); break;
    }
        
    nev=ana->tree()->GetEntries();
    cout << "Selected events="<<nev<<endl; 
    //cout << "NOT STARTING FORM 1ST EVENT!!!!!!!!" << endl;
    for (int i=0; i<nev; i++) {
        //if (!event( ev[i]) ) continue;
        //if (!event(-1 )) continue;
        event(-1);
        //analyze("y");
        analyze("x"); // recoils
    }
    TFile f("segout.root","recreate");
    ana->nt->Write();
    f.Write();
}



TString cut0="wire0.y>wire1.y&&wire0.widthR-wire0.widthL>5&&wire0.y<1e5&&wire0.widthL>5&&wire0.widthR<91";
TString cut1="wire0.y<wire1.y&&wire1.widthR-wire1.widthL>5&&wire1.y<1e5&&wire1.widthL>5&&wire1.widthR<91";

void plotEnergyRange(TString opt="") {
        ana->nt->Draw("(wire0.widthR-wire0.widthL)/58.8*8:wire0.y/10", cut0, opt);
        ana->nt->Draw("(wire1.widthR-wire1.widthL)/58.8*8:wire1.y/10", cut1, opt+"same");
}

void plotEnergyRangeSRIM(TString what="F long", TString opt="L") {
        opt.ToLower();
        MaxCamMC mc;
        mc.setProjectile(14.1e3, 0, 0, 1e6); // 14.1 MeV neutron (all in keV)
        if (what.Contains("F")) { mc.fillSrimTable("SRIM_F_in_CF4_100Torr"); mc.setRecoil(0,0,0,19e6); } // F -> CF4  
        if (what.Contains("C")) { mc.fillSrimTable("SRIM_C_in_CF4_100Torr"); mc.setRecoil(0,0,0,12e6); } // C -> CF4  
        mc.setStopping(200); // Torr
        if      (what.Contains("long")) mc.getRangeVsEnergyProject("long")->Draw(opt);
        else if (what.Contains("tran")) mc.getRangeVsEnergyProject("tran")->Draw(opt);
        else  mc.getRangeVsEnergy()->Draw(opt);
}

void plotEnergyRangeAlphaSRIM(TString opt="L") {
        MaxCamMC mc;
        mc.fillSrimTable("SRIM_He_in_CF4_100Torr"); // He -> CF4
        mc.setStopping(200); //  Torr
        mc.getRangeVsEnergy()->Draw(opt);
}

void plotSkewness() {
        ana->nt->Draw("wire0.skewness:(wire0.widthR-wire0.widthL)", cut0);
        ana->nt->Draw("wire1.skewness:(wire1.widthR-wire1.widthL)", cut1, "same");
}

