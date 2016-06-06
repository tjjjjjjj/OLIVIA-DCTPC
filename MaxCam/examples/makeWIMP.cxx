#include "../MaxCamMC.hh"
#include "../MaxCamWIMP.hh"
#include "TF1.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"

#include <iostream>
using std::cout;
using std::endl;

double wimpMass=100e6; // in keV

MaxCamMC *mc=0;


int nbinning=-1;

void begin(TString fname="mcrun.root") {
    gSystem->Setenv("MCTABLES","../tables");

    
    mc=new MaxCamMC(fname);
    
    int n=6;
    mc->setTrackImage(96*12/n,0,46.7,  64*12/n,0,31.2); // chamber roi in mm (same as 'track image')
    mc->setAnodeImage(96*12/n,0,46.7,  64*12/n,0,31.2); // chamber roi in mm (same as 'track image')
    mc->setCCDImage  (96*12/n,0,768,   64*12/n,0,512); // ccd coordinates (same binning as 'wire image')
    mc->fillSrimTable("SRIM_F_in_CF4_100Torr"); // SRIM table for F -> CF4
    mc->setPressure(75);  // Torr
    mc->setPhotonsPerkeV( 100 /*27*/ ); //  counts/keV
    mc->setPixelsPermm(96/46.7 *12/n); // bins/mm
    mc->setRecoil(0,0,0,19e6); // set recoil nucleus to fluorine
    mc->setNuclScintillation(0.3); // scintillation based on nuclear stopping power
    mc->setNoiseADC( nbinning /*25*/); // 25
   
    mc->clearWireList();
}


void event(float e, bool forceEvent=false) {    
    double MD=wimpMass;

    while (1) {
        // first generate WIMP velocities
        TVector3 v(244,0,0);
        MaxCamWIMP::generateMaxwellWIMP(v);
        double   f=0.5*MD/pow(3e5,2)*v.Mag2();
        TVector3 p=v.Unit();
        p *= f;
        mc->setProjectile(p.X(), p.Y(), p.Z(), wimpMass);
        
        // set recoil energy, check if kinematics possible...
        //mc->setRecoilCoord(25,15,0);// fixed interaction point
        mc->setRandomRecoilCoord(25,35, 15,20, 0,250);
        if (!mc->setRecoilEnergy( e )) {
            // then set vertex point and propagate track
            mc->event(); // make event, save to ntuple
            break;
        }
        
        if (!forceEvent) break; 
    }
    
    
    // display
    //mc->getTrackImage()->Draw("colz");
    //mc->getAnodeImage()->Draw("colz");
    mc->getCCDImage()->Draw("colz");
    
}


void end() {
    mc->Write();
}


void makeWIMP(TString fname="mcrun.root") {
    begin(fname);
    //for (int e=20; e<201; e+=20) {
    //    cout << "Energy: " << e << endl;
    //    for (int i=0; i<1; i++) event(e); //25000
    //}

    // uniform distribution in Erecoil
    //float e;
    //for (int i=0; i<5000; i++) { // 40000 
    //    if (i%100==0) cout << mc->tree()->GetEntries() << "  " << flush;
    //    e=20+gRandom->Rndm()*180;
    //    event(e, true); 
    //}

    // erecoil=100keV
    float e;
    for (int i=0; i<2000; i++) { // 40000 
        if (i%100==0) cout << mc->tree()->GetEntries() << "  " << flush;
        e=100;
        event(e, true); 
    }
    
    end();
}






TF1 *fWimp;

double fERecoil(double *x, double *par) {
    if (par);
    return MaxCamWIMP::dRdE(x[0], wimpMass, 19e6)/MaxCamWIMP::dRdE(0, wimpMass, 19e6);
}


TTree *nt=0;
TProfile *h;
TH1F *hf;
void tmp() {
    begin();
    fWimp=new TF1("fWimp",&fERecoil,0, 500, 0);
    double MD=wimpMass;

    gRandom->SetSeed(1);
    TFile *file= new TFile("mcrun.root","recreate");
    nt=new TTree("nt","");
    float var[8];
    nt->Branch("res", &var, "vx/F:vy:vz:px:py:pz:ERecoil:cosRecoil");
    for (int e=10; e<150; e+=10) {
        cout << "Generating energy " << e << endl;
        for (int i=0; i<10000; i++) {
            TVector3 v(244,0,0);
            MaxCamWIMP::generateMaxwellWIMP(v);
            
            double   f=0.5*MD/pow(3e5,2)*v.Mag2();
            TVector3 p=v.Unit();
            p *= f;
            mc->setProjectile(p.X(), p.Y(), p.Z(), MD);
            
            double cosRecoil=-2, ERecoil=-2;
            if (!mc->setRecoilEnergyAngle( e, 1.57)) {
                ERecoil =  mc->getRecoil()->Vect().Mag();
                cosRecoil = mc->getRecoil()->Vect().X()/ERecoil;
            };        

            var[0]=v.X(); var[1]=v.Y(); var[2]=v.z();
            var[3]=p.X(); var[4]=p.Y(); var[5]=p.Z();
            var[6]=ERecoil; var[7]=cosRecoil;
            nt->Fill();            
        }        
    }
    nt->Write();
    file->Write();
    h=new TProfile("h","",14,10,150,"S");
    h->SetFillStyle(4100); h->SetFillColor(4); h->SetMarkerStyle(1); h->SetMinimum(0);
    h->SetYTitle("cosine of recoil angle"); h->SetXTitle("Recoil energy (keV)");
    nt->Draw("abs(res.cosRecoil):res.ERecoil>>h");
    hf=new TH1F("hf","",14,10,150);
    hf->SetLineStyle(2);
    for (float e=10; e<150; e+=10) { hf->Fill(e, fWimp->Eval(e) ); }
}

