// $Id: McDarkEventAction.cc,v 1.30 2011/03/24 18:33:08 jplopez Exp $
// GEANT4 tag $Name:  $
//
 
#include "McDarkEventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4PrimaryParticle.hh"

//#include "G4UImanager.hh"
#include "G4VVisManager.hh"

#include "TH2.h"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"

#include "G4TransportationManager.hh"

#include "McDarkDigiManager.hh"
#include "McDarkTpcDigitizer.hh"
#include "McDarkRunAction.hh"
#include "McDarkTpcHit.hh"
#include "McDarkCameraManager.hh"

#include "McDarkTrack.hh"
#include "McDarkDigi.hh"
#include "DmtpcEvent.hh"
#include "DmtpcDataset.hh"



#include <iomanip>
using std::setw;

McDarkEventAction::McDarkEventAction(McDarkRunAction *run) : _run(run) { }

McDarkEventAction::~McDarkEventAction() { }

void
McDarkEventAction::BeginOfEventAction(const G4Event* evt) {        
    getRun()->data()->event()->setEventNumber( evt->GetEventID() );
}




void
McDarkEventAction::EndOfEventAction(const G4Event* evt) {    


    // Loop over all registered digitization modules and make digis
    McDarkDigiManager::Digitize();

    // export MC to ROOT format
    if (isTriggered(evt)) {
      exportTracks(evt);
      exportCcdInfo(evt);
      getRun()->data()->fill();
      getRun()->data()->clearEventMemory();    
    }
    
    // add event to display
    if (getRun()->getPlotEvent()) evt->Draw();
    if (getRun()->getPrintTracks()>0) printTrackTable(evt);
}



void
McDarkEventAction::printTrackTable(const G4Event* evt) {

    if (getRun()->getPrintTracks()<1) return;
        
    G4HCofThisEvent* allHits=evt->GetHCofThisEvent();
    if (!allHits) return;
    G4DCofThisEvent* allDigis=evt->GetDCofThisEvent();
    if (!allDigis) return;
    G4int nHC=allHits->GetCapacity();
    G4int nDC=allDigis->GetCapacity();        

    
    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    if (!trajectoryContainer) return;


    int n_trajectories = trajectoryContainer->entries();
    if (getRun()->getPrintTracks()==1) { // option '1'
        if (n_trajectories==1) return; // don't print if no interaction
        if (n_trajectories>10) n_trajectories=10; // limit printout to 10 tracks
    }

    G4cout << G4endl;
    G4cout << "Event number = " << evt->GetEventID() << G4endl;
    G4cout << "Index  |         Particle | Mother |     E (keV) |  Range (cm) |  Eion (keV) |   CCD (ADU) |" << endl;
    G4cout << "-------+------------------+--------+-------------+-------------+-------------+-------------|" << endl;

    for ( G4int i=0; i<n_trajectories; i++) {
        
        // Note: entries in track table are not ordered;
        // order them here for easier handling
        G4Trajectory* t=0;
        for (int ii=0; ii<trajectoryContainer->entries(); ii++) {
            t = (G4Trajectory*)((*(trajectoryContainer))[ii]);
            if (t->GetTrackID()==i+1) break;
        }

        // get kin. energy
        G4double Ekin = t->GetInitialKineticEnergy();
        
        // find range
        const G4ThreeVector &begPoint=t->GetPoint(0)->GetPosition();
        const G4ThreeVector &endPoint=t->GetPoint( t->GetPointEntries() - 1 )->GetPosition();
        G4double range=(begPoint-endPoint).mag();

        
        // find quenched energy
        G4double quenchedEnergy=0;
        for (int iHC=0; iHC<nHC; iHC++) {
            McDarkTpcHitCollection *hc = (McDarkTpcHitCollection*)allHits->GetHC( iHC );
            for (unsigned int iHit=0; iHit<hc->GetSize(); iHit++) {
                McDarkTpcHit *hit=(McDarkTpcHit*)hc->GetHit(iHit);
                if (hit->getTrackID()!=t->GetTrackID()) continue;
                quenchedEnergy+=hit->getDepositEnergy()*hit->getQuenchingFactor();
            }
            
        }

        
        // find CCD yield
        G4double ccdYield=0;
        for (int iDC=0; iDC<nDC; iDC++) {
             McDarkTpcDigiCollection *dc = (McDarkTpcDigiCollection*)allDigis->GetDC( iDC );
             G4String dcName = dc->GetName();
             if (dcName.contains("ccd")) {
                for (unsigned int iDigi=0; iDigi<dc->GetSize(); iDigi++) {
                    McDarkTpcDigi *digi=(McDarkTpcDigi*)dc->GetDigi(iDigi);
                    if (digi->getTrackID()!=t->GetTrackID()) continue;
                    ccdYield+=digi->getWeight();
                }
            }
        }
            

        
        // print track info
        G4cout.precision(5);
        G4cout << setw(6)  << t->GetTrackID() << " | "
               << setw(16) << t->GetParticleName() << " | "                   
               << setw(6) << t->GetParentID() << " | "
               << setw(11) << Ekin/keV << " | "
               << setw(11) << range/cm << " | "
               << setw(11) << quenchedEnergy/keV << " | "
               << setw(11) << ccdYield << " | "
               << G4endl; 
    }
    

    if ( n_trajectories < trajectoryContainer->entries() ) {
        G4cout << "(Shown " << n_trajectories
               << " of " << trajectoryContainer->entries()
               << " tracks)"
               << G4endl;
    }


}


G4bool
McDarkEventAction::isTriggered(const G4Event* evt) {

    // require digitized hits in detector
    G4DCofThisEvent* allDigis=evt->GetDCofThisEvent();
    if (!allDigis) return false;
    G4int nDC=allDigis->GetCapacity();        
    for (int iDC=0; iDC<nDC; iDC++) {
        McDarkTpcDigiCollection *dc = (McDarkTpcDigiCollection*)allDigis->GetDC( iDC );
        if (dc && dc->GetSize()>0) return true;
    } 
    return false;
}

void
McDarkEventAction::exportTracks(const G4Event* evt) {

    // get number of stored trajectories
    G4int n_trajectories = 0;
    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

    // access to hits
    G4HCofThisEvent* allHits=evt->GetHCofThisEvent();
    G4int nHC=allHits->GetCapacity();

    
    int iSaved=0;
    for ( G4int i=0; i<n_trajectories && i<DmtpcEvent::MAX_N_MCDARK_TRACKS; i++) {

        // Note: entries in track table are not ordered;
        // order them here for easier handling
        G4Trajectory* t=0;
        for (int ii=0; ii<n_trajectories; ii++) {
            t = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[ii]);
            if (t->GetTrackID()==i+1) break;
        }
        
        const G4ThreeVector &begPoint=t->GetPoint(0)->GetPosition();
        const G4ThreeVector &endPoint=t->GetPoint( t->GetPointEntries() - 1 )->GetPosition();     

        // find quenched energy
        G4double quenchedEnergy=0;
        for (int iHC=0; iHC<nHC; iHC++) {
            McDarkTpcHitCollection *hc = (McDarkTpcHitCollection*)allHits->GetHC( iHC );
            for (unsigned int iHit=0; iHit<hc->GetSize(); iHit++) {
                McDarkTpcHit *hit=(McDarkTpcHit*)hc->GetHit(iHit);
                if (hit->getTrackID()!=t->GetTrackID()) continue;
                quenchedEnergy+=hit->getDepositEnergy()*hit->getQuenchingFactor();
            }            
        }
        
        const G4ThreeVector &p=t->GetInitialMomentum();
        double m=t->GetParticleDefinition()->GetPDGMass();
        McDarkTrack aTrack;
        aTrack.setP4( p.x()/keV, p.y()/keV, p.z()/keV, sqrt( p.mag2() + m*m )/keV );
        aTrack.setQuenchedEnergy( quenchedEnergy/keV );
        aTrack.setIndex(t->GetTrackID());
        aTrack.setMotherIndex(t->GetParentID());
        aTrack.setBegPoint( TVector3( begPoint.x(), begPoint.y(), begPoint.z() ) );
        aTrack.setEndPoint( TVector3( endPoint.x(), endPoint.y(), endPoint.z() ) );
        aTrack.setPDGIndex( t->GetParticleDefinition()->GetPDGEncoding() );
        new( (*getRun()->data()->event()->mcTrack())[iSaved++] ) McDarkTrack( aTrack );
    }    
        
}




void
McDarkEventAction::exportCcdInfo(const G4Event *evt) {

    // Save all CCD digis saved in digi collection of event
    G4DCofThisEvent* allDigis=evt->GetDCofThisEvent();
    G4int nDC=allDigis->GetCapacity();
    G4int iSaved=0;
    for (int iDC=0; iDC<nDC; iDC++) {
        McDarkTpcDigiCollection *dc = (McDarkTpcDigiCollection*)allDigis->GetDC( iDC );
        G4String dcName = dc->GetName();
        if (dcName.contains("ccd")) {
            for (unsigned int iDigi=0; iDigi<dc->GetSize(); iDigi++) {
                McDarkTpcDigi *digi=(McDarkTpcDigi*)dc->GetDigi(iDigi);
                McDarkDigi aDigi;
                aDigi.setModuleID( digi->getModuleID() );
                aDigi.setChannelID( digi->getChannelID() );
                aDigi.setWeight( digi->getWeight() );
                aDigi.setTrackID( digi->getTrackID() );
                new( (*getRun()->data()->event()->mcCcdDigi())[iSaved++] ) McDarkDigi( aDigi );
            }
        }
    }
    

    // Save images from digitization modules
    iSaved=0;
    unsigned int ncam = McDarkCameraManager::listOfCameras.size();
    for (unsigned int icam=0; icam<ncam; icam++) { 
    
      McDarkCamera *cam=McDarkCameraManager::listOfCameras[icam];
      cam->addNoise();
      TH2S* im = cam->rawImage();
      new( (*getRun()->data()->event()->ccdData())[iSaved] ) TH2S( *im  );
      new( (*getRun()->data()->event()->ccdConfig())[iSaved] ) MaxCamConfig( *cam->config()  );
      cam->clearImage();
      iSaved++;
      delete im;
    }

    
}
