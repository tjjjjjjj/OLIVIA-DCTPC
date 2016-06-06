//--*-C++-*--

#include "../MaxCamMC.hh"
#include "TF1.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"

#include <iostream>
using std::cout;
using std::endl;

MaxCamMC *mc=0;

void begin() {
    gSystem->Setenv("MCTABLES","/Users/ddujmic/work/projects/DarkMatter/MaxCam/tables");
    gRandom->SetSeed(3000);

    mc=new MaxCamMC("mcrun.root"); 

    mc->setWireImage(96,0,12,  64,-1.4,6.6); // chamber roi in mm x mm
    mc->setCCDImage (96,0,768,  64,0,512); // ccd coordinates (same binning as 'wire image')
    mc->fillSrimTable("SRIM_F_in_CF4_100Torr"); // SRIM table for F -> CF4
    mc->setPressure(200);  // Torr
    mc->setPhotonsPerkeV(8); //  counts/keV
    mc->setPixelsPermm(96/12.); // bins/mm
    mc->setProjectile(-14.1e3, 0, 0, 1e6); // 14.1MeV neutron
    mc->setNuclScintillation(0.3); // generate scintillation based on elec stopping power
    mc->setNoiseADC(25);

    mc->clearWireList();
    for (int i=0; i<20; i++) mc->addWire( -50+i*5 ); // wires in mm
}




int event(float e) {
        
    // create event
    //mc->setRecoilEnergy(e); // set recoil energy in keV
    // or
    //mc->setRecoilEnergyAngle(e,1.57); // recoil energy + phi angle
    // or
    mc->setRecoil(-e, 0, 0, 19e6); // set a track with these parameters

    //return ( mc->setRecoilEnergy(e) );

    
    //mc->setRandomRecoilCoord(-4,16, -5.4,9); // interaction volume/plane
    mc->setRecoilCoord(10, 0, 20);
    mc->event(); // make event
    
    // display
    mc->getTrackImage()->Draw("colz");
    //mc->getWireImage()->Draw("colz");
    //mc->getCCDImage()->Draw("colz");

    return 0;
}


void end() {
    mc->Write();
}


void all() {
    begin();
    for (float e=25; e<150; e+=25)
        for (int i=0; i<1000; i++) event(e);
    end();
}





TNtuple *nt, *ntsrim;
void tmp() {
    begin();
    nt= new TNtuple("nt","","x:y:z");
    for (int i=0; i<1000; i++) {
        event(1000);
        nt->Fill( mc->getRecoilCoord()->X(), mc->getRecoilCoord()->Y(), mc->getRecoilCoord()->Z()  );
    }
}


void tmpsrim() {
    ifstream fin("RANGE_3D.txt");    
    string line;
    float x,y,z;
    ntsrim= new TNtuple("nt","","x:y:z");
    while (!fin.eof()) {
        getline( fin, line);
        istringstream sline(line);
    	sline >> event >> x >> y >> z;
        ntsrim->Fill(x*1e-7,y*1e-7,z*1e-7);
    }
}


void energyangle() {
    begin();
    TFile *file= new TFile("mcrun.root","recreate");
    nt= new TNtuple("nt","","px:py:pz:e");
    for (float e=5; e<150; e+=5) {
        for (int i=0; i<100; i++) {
            if (event(e)) continue;
            nt->Fill( mc->getRecoil()->X(), mc->getRecoil()->Y(), mc->getRecoil()->Z(), e  );
        }
    }
    end();
    nt->Write();
    file->Write();
}
