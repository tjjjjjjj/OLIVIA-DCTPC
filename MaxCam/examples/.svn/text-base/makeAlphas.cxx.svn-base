#include "../MaxCamMC.hh"
#include "../MaxCamENDF.hh"
#include "../MaxCamWIMP.hh"
#include "../MaxCamTwoBodyKinematics.hh"
#include "TF1.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TStyle.h"

#include <iostream>
using std::cout;
using std::endl;


MaxCamMC *mc=0;

void begin() {
    gSystem->Setenv("MCTABLES","/Users/ddujmic/work/projects/DarkMatter/MaxCam/tables");
    
    mc=new MaxCamMC("mcrun.root"); 

    mc->setWireImage(96,0,46.7,  64,0,31.2); // chamber roi in mm (same as 'track image')
    mc->setCCDImage (96,0,768,   64,0,512); // ccd coordinates (same binning as 'wire image')
    mc->fillSrimTable("SRIM_He_in_CF4_100Torr"); // SRIM table for He -> CF4
    mc->setPressure(75);  // Torr
    mc->setPhotonsPerkeV(27); //  counts/keV
    mc->setPixelsPermm(96/46.8); // bins/mm
    mc->setProjectile(-500, 0, 0, 4e6); // alpha direction, mass
    mc->setNuclScintillation(0.3); // scintillation based on nuclear stopping power
    mc->setNoiseADC(25);
        
    mc->clearWireList();
    //for (float d=-50; d<50; d+=5) mc->addWire(d);

}



 



void event(double e=2000) {
        
    // create event

    
    double phi=0;//TMath::Pi()/4;
    double theta=2*TMath::Pi()*gRandom->Rndm();
    mc->setRecoil(-e*cos(phi), e*sin(theta)*sin(phi), e*cos(theta)*sin(phi), 4e6);

    mc->setRecoilCoord(30,8,0);// fixed interaction point

    mc->event(); // make event

    // display
    //mc->getWireImage()->Draw("colz");
    //mc->getTrackImage()->Draw("colz");
    mc->getCCDImage()->Draw("colz");

}




void end() {
        mc->Write();
}


void all() {
    begin();
    for (float e=20; e<940; e+=20) {
        cout << "Generate energy: " << e << "keV"<<endl;
        for (int i=0; i<1000; i++) event(e);
    }
    end();
}

