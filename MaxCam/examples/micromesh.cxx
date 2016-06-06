//////////////////////////////////////////////////////////////////////////////////
//
//   Recoinstruction of recoil segments for neutron scattering (DT) experiment
//   Data runs 139+146+149
//   Publications: NIM/EPS conf
//
//   - in ROOT to run on run number NNN: 
//   .L segments.cxx++
//   segments(NNN)
//   - output file is RESULTS/results_runNNN.root
//   - make plots using segmentsPlot_DT.cxx
//
//////////////////////////////////////////////////////////////////////////////////


#include "../MaxCamRead.hh"
#include "../MaxCamMC.hh"
#include "../MaxCamTrack.hh"
#include "../MaxCamImageTools.hh"
#include "../MaxCamSegment.hh"
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


void init(int runN, int outRun=-1) {
  //int col[]={10, 19,18,17,16,15,14,13,12, 1};
  //gStyle->SetPalette(sizeof(col)/sizeof(int),col);
  
  // OLD   gSystem->Setenv("MCTABLES","/Users/ddujmic/work/projects/DarkMatter/MaxCam/tables");
  gSystem->Setenv("MCTABLES","../tables");

  TString name="data/ccdrun_00"; name+=runN; name+=".root";
  cout << " *************** Input "<< name << " ***************** " << endl;
  TString outputFileName="";
  if (outRun>0) {
      outputFileName="RESULTS/results_run"; outputFileName+=outRun; outputFileName+=".root";
      cout << " *************** Output "<< outputFileName << " ***************** " << endl;
  }
  ana = new MaxCamRead(name, outputFileName);   
  ana->readBiasFrame();
  ana->findHotPixels("hotpixels.dat");
  
  //TString name="MonteCarloSamples/mcrun_000"; name+=runN; name+=".root";
  //cout << " *************** Input "<< name << " ***************** " << endl;
  //TString outputFileName="RESULTS/results_mcrun_000"; outputFileName+=runN;  outputFileName+=".root";
  //cout << " *************** Output "<< outputFileName << " ***************** " << endl;
  //ana = new MaxCamRead(name, outputFileName);
  

  //ana = new MaxCamRead("mcrun.root", "results_mcrun.root");

  
  c1=new TCanvas("c1","",0,0,700, 600); 
  c2=new TCanvas("c2","",700,0,750, 600);
  c1->cd();
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
}

#include "TList.h"
#include "TPolyMarker.h"



int analyze(TString what) {
    what.ToLower();

    
    ana->segment()->reset();
    if (what=="cf") { // neutron scattering
        c2->cd();
        
        //////////////////////////////////////////////////////////////////////////////////////////////if (ana->isDischargeEvent(100)) return 0;
        
        int imax, jmax;
        float threshold;
        if (!ana->hasSegments(imax, jmax, threshold)) {
            if (ana->isMC()) fill(); // no recoil?  save mc truth only
            return 0; // for data do not save event
        }
        ana->ccdImage()->DrawCopy("colz");

        TH2F* cleanImage=(TH2F*)ana->ccdImage()->Clone("cleanImage");

        
        c1->cd();
        ana->getEvent(i-1);
        ana->ccdImage()->Draw("colz");
        float bg=0;
        double I = MaxCamImageTools::calcIntensity2D(ana->ccdImage(), imax, jmax, 5, bg);
        ana->segment()->setEnergy( I, 0 );

        ana->segment()->setSkewness( MaxCamImageTools::calcSkewness2D(cleanImage, imax, jmax, bg) );
        delete cleanImage;
        
        ana->segment()->print();
    }

    else if (what=="am") { // alpha calibration
        
        // tracking
        TH1F *hY = ana->makeYieldHisto();
        MaxCamTrack* trfit = new MaxCamTrack( ana->ccdImage(), true );
        float threshold = hY->GetMean() + hY->GetRMS()*3;
        trfit->setThreshold( threshold );
        trfit->makeTracks(); 
        ana->segment()->setNSegments( trfit->nTracks() );
        //cout << ": Found # of tracks="<< trfit->nTracks() << endl;
        //if (trfit->nTracks()) return false; // recoils
        float slope = trfit->nTracks() ? trfit->getTrack(0)->GetParameter(1) : 0;
        //if (fabs(slope)>0.1) return false;
        //cout << "slope="<<trfit->getTrack(0)->GetParameter(1)<< endl;
        ana->segment()->setCorrelation( slope );
        delete hY;
        delete trfit;    
        
        // count 
        c2->cd();
        TH1D * hy0 = ana->ccdImage()->ProjectionY("_py", 2, 30);
        //hy0->GetXaxis()->SetRange( 1, 29); // upper drift
        //hy0->GetXaxis()->SetRange(36, 96); // lower drift
        //TH1D * hy0 = ana->ccdImage()->ProjectionX("_py", 50, 60);
        hy0->ShowPeaks(1, "", 0.5);
        TList *functions = hy0->GetListOfFunctions();
        TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
        int nseg=0;
        for (int jj=0; pm && jj<pm->GetN(); jj++) {
            if (pm->GetY()[jj]>1000) nseg++;
            //cout << pm->GetX()[jj]<<"  " << pm->GetY()[jj] << endl;
        }
        ana->segment()->setNSegments( nseg );

        // energy
        ana->segment()->setDischarge( ana->isDischargeEvent(100) ? true : false );
        TH1D * hy = ana->ccdImage()->ProjectionY("_py", 2,12);
        hy->GetXaxis()->SetRange( 5, 56); // upper drift
        //hy->GetXaxis()->SetRange(36, 96); // lower drift
        //TH1D * hy = ana->ccdImage()->ProjectionX("_py", 75,85);
        hy->Draw();
        ana->getProfileFunction()->SetParameters(hy->GetMinimum(), hy->GetMaximum(), hy->GetBinCenter(hy->GetMaximumBin()), 3);
        hy->Fit("_fun");
        TF1* tfun=hy->GetFunction("_fun");
        double I = sqrt(2*3.14)*tfun->GetParameter(1)*tfun->GetParameter(3)/hy->GetBinWidth(1);
        double IErr = I * sqrt( pow( tfun->GetParError(1)/tfun->GetParameter(1),2) +
                                pow( tfun->GetParError(3)/tfun->GetParameter(3),2 ) );
        ana->segment()->setEnergy(I, IErr);
        ana->segment()->setBackground(  tfun->GetParameter(0), tfun->GetParError(0) );
        ana->segment()->setAmplitude(   tfun->GetParameter(1), tfun->GetParError(1) );
        ana->segment()->setMean(        tfun->GetParameter(2), tfun->GetParError(2) );
        ana->segment()->setWidth(       tfun->GetParameter(3), tfun->GetParError(3) );
    }
    else if( what == "fe" ) {
        if (ana->isDischargeEvent(100)) return 0;
        c2->cd();
        for (int ii=1; ii<=93; ii++) {
            for (int jj=26; jj<=26; jj++) {
                double neighbor=0;
                for (int k=ii-1; k<=ii+1; k++) for (int l=jj-1;l<=jj+1; l++) neighbor+=ana->ccdImage()->GetBinContent(k,l);
                ana->segment()->setEnergy( ana->ccdImage()->GetBinContent(ii,jj), neighbor ); 
                ana->segment()->setNPixels(ii);
                fill(); 
            }
        }
        return 0;
    }
    
    fill();
    return 0;
}



void segments(int runN, TString what) {

    // what =  am, fe, cf
    

    if (!ana) init(runN, runN);
    what.ToLower();
    
    i=0;
    int nev=ana->tree()->GetEntries();
    cout << "Selected events="<<nev<<endl; 
    //cout << "NOT STARTING FORM 1ST EVENT!!!!!!!!" << endl;
    for (; i<nev; i++) {
        //if (!event(-1 )) continue;
        event(-1);
        analyze(what); 
    }

    ana->closeOutputFile();
}



