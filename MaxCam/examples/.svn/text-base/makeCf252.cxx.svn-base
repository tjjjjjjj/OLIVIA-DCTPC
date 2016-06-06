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
MaxCamENDF *neutronEnergy=0, *scatteringDCS=0, *scatteringCS=0;

//TString ionType="4He";
TString ionType="19F";
double ionMass=0, projMass=1e6;

void begin() {
    gSystem->Setenv("MCTABLES","/Users/ddujmic/work/projects/DarkMatter/MaxCam/tables");


    neutronEnergy=new MaxCamENDF("ENDF_Cf-252_n_spectrum", "fission");
    scatteringDCS=new MaxCamENDF(TString("ENDF_DCS_n_on_")+ionType, "elastic scattering");
    scatteringCS = new MaxCamENDF(TString("ENDF_CS_n_on_")+ionType, "cs");    
    
    mc=new MaxCamMC("mcrun.root"); 
    
    mc->setRecoil(0,0,0,4e6); 
    if (ionType=="19F") {
        mc->fillSrimTable("SRIM_F_in_CF4_100Torr"); // SRIM table for F -> CF4
        ionMass=19e6; //keV
    }
    else if (ionType=="4He") {
        mc->fillSrimTable("SRIM_He_in_CF4_100Torr"); // helium
        ionMass=4e6; //keV
    }
    mc->setRecoil(0,0,0,ionMass);
    projMass=1e6;

    // set up detector:
    mc->setWireImage(96,0,46.7,  64,0,31.2); // chamber roi in mm (same as 'track image')
    mc->setCCDImage (96,0,768,   64,0,512);  // ccd coordinates (same binning as 'wire image')
    mc->setPressure(75);  // Torr
    mc->setPhotonsPerkeV(27); //  counts/keV
    mc->setPixelsPermm(96/46.8); // bins/mm
    mc->setNuclScintillation(0.3); // scintillation based on nuclear stopping power
    mc->setNoiseADC(25);
    
    mc->clearWireList(); // mesh detector
     
}

void event(float emin, float emax) {

    double e=0;
    while (1) {
        double Eneutron=neutronEnergy->generateEnergy(); // neutron energy distribution from Cf-252
        mc->setProjectile(-Eneutron, 0, 0, projMass);

        if (!scatteringCS->acceptEnergy(Eneutron)) continue; // cross section
        
        double cosScatterCMS=scatteringDCS->generateCosAngleCMS( Eneutron );
        e = MaxCamTwoBodyKinematics::calcRecoilEnergyFromCosScatterCMS(cosScatterCMS, Eneutron, ionMass, projMass);

        if (e>emin && e<emax) break;
    }
    cout << "energy = " << e << endl;
    
    
    mc->setRecoilCoord(25,15,0);// fixed interaction point
    
    if ( !mc->setRecoilEnergy(e) ) {  
        mc->event(); // make event
     }
    
    
    // display
    //mc->getTrackImage()->Draw("colz");
    //mc->getAnodeImage()->Draw("colz");
    //mc->getCCDImage()->Draw("colz");
}




void event(float e) {

    double Eneutron=neutronEnergy->generateEnergy(); // neutron energy distribution from Cf-252
    mc->setProjectile(-Eneutron, 0, 0, 1e6);


    
    mc->setRecoilCoord(25,15,0);// fixed interaction point
    
    if ( !mc->setRecoilEnergy(e) ) {  
        mc->event(); // make event
     }
    
    
    // display
    //mc->getTrackImage()->Draw("colz");
    //mc->getAnodeImage()->Draw("colz");
    mc->getCCDImage()->Draw("colz");
}


void end() {
    mc->Write();
}


void makeCf252() {
    begin();
    //for (int e=20; e<941; e+=20) for (int i=0; i<500; i++) event(e);
    for (int i=0; i<10000; i++) event(10, 1000);
    end();
}


TNtuple *ntp;
void espec() {
    begin();

    ntp=new TNtuple("ntp","","en:cosRecoil:er");
    for (int i=0; i<10000; i++) {
        double Eneutron=neutronEnergy->generateEnergy(); // neutron energy distribution from Cf-252
        if (!scatteringCS->acceptEnergy(Eneutron)) continue; // cross section
        //cout << Eneutron << endl;
        double cosScatter=scatteringDCS->generateCosAngleCMS(Eneutron); // recoil energy distribution
        //cout << cosScatter << endl;
        double Erecoil= MaxCamTwoBodyKinematics::calcRecoilEnergyFromCosScatterCMS( cosScatter, Eneutron, ionMass, projMass);
        double cosRecoil = MaxCamTwoBodyKinematics::calcCosRecoilFromRecoilEnergy(Erecoil, Eneutron, ionMass, projMass) ;

        ntp->Fill(Eneutron,cosRecoil, Erecoil);
        
    }
    end();
}

