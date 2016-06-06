//
//-------------------------------------------------------------------
//
// Created by: Timur Sahin (tcsahin@MIT.EDU)
// Date:       July 15, 2007
// Copyright:  MIT 2007
//
//
//
//
// $Id: McDarkPrimaryGeneratorAction.cc,v 1.27 2011/01/07 15:42:58 ddujmic Exp $
// GEANT4 tag $Name:  $
//

#include "McDarkPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <assert.h>

#ifndef G4NOHIST
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#endif

#include "McDarkSpergelDistribution.hh"




McDarkPrimaryGeneratorAction::McDarkPrimaryGeneratorAction() {
    
    particleGun = new G4ParticleGun();

    spergelDist = new McDarkSpergelDistribution();

    theMessenger = new McDarkPrimaryGeneratorActionMessenger(this);

 
}



McDarkPrimaryGeneratorAction::~McDarkPrimaryGeneratorAction()
{
  delete particleGun;
  delete spergelDist;
}


void McDarkPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

    const char* seedString = gSystem->Getenv("MCDARKSEED");
    if (seedString) {
      int theseed = atoi(seedString);
      CLHEP::HepRandom::setTheSeed(theseed);
    }

    if ( distributionType==Spergel ) { // Spergel distribution        
      G4cout << "This is not spergeljob" << G4endl;

        int A = particleGun->GetParticleDefinition()->GetAtomicNumber();
        getSpergelFunction()->setTargetAtomicNumber(A);
        
        Double_t e_and_psi[2];
        spergelDist->GetRandomValues(e_and_psi);
        G4double en = e_and_psi[0];
        G4double psi = e_and_psi[1];
    
        G4ThreeVector v(TMath::Cos(psi), TMath::Sin(psi), 0);
        Double_t phi = G4UniformRand()*360;
        v.rotateX(phi*deg);

        particleGun->SetParticleEnergy(en*MeV);
        particleGun->SetParticleMomentumDirection(v);
	particleGun->GeneratePrimaryVertex(anEvent);
    }

    else if ( distributionType==Gun ) { //  simple particle gun with predefined parameters

        G4double xRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minX, particle_maxX);
        G4double yRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minY, particle_maxY);
        G4double zRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minZ, particle_maxZ);
        G4ThreeVector dir(xRand, yRand, zRand);
        particleGun->SetParticlePosition(dir);
       
        double energy = G4UniformRand()*(particle_maxEnergy-particle_minEnergy)+particle_minEnergy; 
        particleGun->SetParticleEnergy( energy );  
        
        particleGun->GeneratePrimaryVertex(anEvent);

    } else if (distributionType == Cf252) {

        // Cf-252 neutron spectrum
        G4double energy = 0, rand=0, fun=0;
        while (1) {
            energy = G4UniformRand()*(particle_maxEnergy-particle_minEnergy)+particle_minEnergy; 
            rand = G4UniformRand();
            fun = exp(-0.93*energy)*sinh(sqrt(1.25*energy));
            //G4cout << GetName() << energy << "  " << rand << "  " << fun << G4endl;
            if ( fun>rand ) break;
        }
        particleGun->SetParticleEnergy( energy );  
       
        G4double xRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minX, particle_maxX);
        G4double yRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minY, particle_maxY);
        G4double zRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minZ, particle_maxZ);
        G4ThreeVector dir(xRand, yRand, zRand);
        particleGun->SetParticlePosition(dir);
        
        particleGun->GeneratePrimaryVertex(anEvent);

    } else if (distributionType == Isotropic) {
        

	G4double rndm, rndm2;
        G4double px, py, pz;

        G4double cosMinTheta=-1, cosMaxTheta=1,  MinPhi=0, MaxPhi=twopi;
        G4double sintheta, sinphi, costheta, cosphi;
        rndm = G4UniformRand();
        costheta = cosMinTheta + rndm * (cosMaxTheta - cosMinTheta);
        sintheta = std::sqrt(1. - costheta*costheta);
        
        rndm2 = G4UniformRand();
        G4double Phi = MinPhi + (MaxPhi - MinPhi) * rndm2; 
        sinphi = std::sin(Phi);
        cosphi = std::cos(Phi);

        px = sintheta * cosphi;
        py = sintheta * sinphi;
        pz = costheta;

        //G4cout <<  G4ThreeVector(px, py, pz) << G4endl;
        particleGun->SetParticleMomentumDirection( G4ThreeVector(px, py, pz) );

	G4double energy = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minEnergy, particle_maxEnergy);

	
        G4double xRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minX, particle_maxX);
        G4double yRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minY, particle_maxY);
        G4double zRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minZ, particle_maxZ);
        G4ThreeVector dir(xRand, yRand, zRand);
        particleGun->SetParticlePosition(dir);
	
	
	
	particleGun->SetParticleEnergy(energy);
	particleGun->GeneratePrimaryVertex(anEvent);
    } 

    else if ( distributionType==DT ||  distributionType==DD) { //  DT or DD soruce

        G4double xRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minX, particle_maxX);
        G4double yRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minY, particle_maxY);
        G4double zRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minZ, particle_maxZ);
        G4ThreeVector dir(xRand, yRand, zRand);
        particleGun->SetParticlePosition(dir);

        if (distributionType==DT) particleGun->SetParticleEnergy(14.1*MeV);
        else  particleGun->SetParticleEnergy(2.6*MeV);
        
        particleGun->GeneratePrimaryVertex(anEvent);

    }

    else if ( distributionType==CosmicMuons ) { // from Bogdanova et al, arXiv:0601019


        double cosMu=0, pMu=0, I=0;
        while (1) {
            cosMu=CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), -1, 1);
            pMu=CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), 1, 1e5);
            I=18./(pMu*cosMu + 145)*pow(pMu+2.7/cosMu, -2.7)*(pMu+5)/(pMu+5/cosMu);
            break;
        }
        assert(0);

    }

    else if (distributionType == Co57) {

      // Co-57 gamma distribution

      //               Fe    14keV 122keV 136keV
      G4double rate[]={0.58, 0.09, 0.11,  0.86};
      G4int nr=sizeof(rate)/sizeof(G4double);
      G4double rateSum=0;
      for (G4int i=0; i<nr; i++) rateSum+=rate[i];
      
      
      G4double rand = G4UniformRand();
      
      
      G4double energy = 0;
      if (rand<rate[0]/rateSum) {
	
	// throw another number and split 50% to 'continuum' and 50% to 6.4 keV peak
	if(G4UniformRand()<0.75) energy =  CLHEP::RandExponential::shoot( 15.0 ) * keV;
	else energy = 6.4 * keV;
	
      } else if (rand < (rate[0]+rate[1])/rateSum) {

	energy = 14.4*keV; 

      } else if (rand < (rate[0] + rate[1] + rate[2])/rateSum) {

	energy = 122.06 * keV;

      } else {

	energy = 136.47 * keV;

      }

      particleGun->SetParticleEnergy( energy );  
      
      G4double xRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minX, particle_maxX);
      G4double yRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minY, particle_maxY);
      G4double zRand = CLHEP::RandFlat::shoot(CLHEP::HepRandom::getTheEngine(), particle_minZ, particle_maxZ);
      G4ThreeVector dir(xRand, yRand, zRand); 
      particleGun->SetParticlePosition(dir);
      
      
      particleGun->GeneratePrimaryVertex(anEvent);

    }
    
    else assert(0);


    print();
}


void
McDarkPrimaryGeneratorAction::print() {

    static G4int isShown=0;

    if (isShown) return;

    G4cout << "+----------------------------------------------------+" << G4endl;
    G4cout << "|                                                    |" << G4endl;
    G4cout << "|  Primary particle generator                        |" << G4endl;
    G4cout << "|                                                    |" << G4endl;
    G4cout << "|       ";
    if ( distributionType==Spergel ) G4cout << GetName() << ": Spergel distribution" << G4endl;
    else if ( distributionType==Gun ) G4cout << GetName() << ": Simple Gun" << G4endl;
    else if (distributionType == Cf252) G4cout << GetName() << ": Cf252 source" << G4endl;
    else if (distributionType == Isotropic) G4cout << GetName() << ": Isotropic source" << G4endl;
    else if (distributionType==DT) G4cout << GetName() << ": DT generator"<<G4endl;
    else if (distributionType==DD) G4cout << GetName() << ": DD generator"<<G4endl;
    else if ( distributionType==CosmicMuons ) G4cout << GetName() << ": Cosmic muons"<<G4endl;
    else if ( distributionType==Co57 ) G4cout << GetName() << ": Co57 source"<<G4endl;
	else assert(0);
    G4cout << "|                                                    |" << G4endl;
    G4cout << "|                                                    |" << G4endl;
    G4cout << "------------------------------------------------------" << G4endl;

    isShown++;
}
