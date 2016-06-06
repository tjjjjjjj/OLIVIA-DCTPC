#include "McDarkTpcDigitizer.hh"
#include "McDarkTpcHit.hh"
#include "McDarkTpcGas.hh"
#include "McDarkCamera.hh"
#include "McDarkCameraManager.hh"

#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4ios.hh"
#include "G4VSensitiveDetector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Tubs.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleTable.hh"

McDarkTpcDigitizer::McDarkTpcDigitizer(G4String name, const G4PVPlacement *det)
    : G4VDigitizerModule(name),
      _detector(det) {

    
    G4String colName = name+"_ccdCollection";
    collectionName.push_back(colName);
    
    colName = name+"_chargeCollection";
    collectionName.push_back(colName);
    
    colName = name+"_pmtCollection";
    collectionName.push_back(colName);  
}


McDarkTpcDigitizer::~McDarkTpcDigitizer() {}


void
McDarkTpcDigitizer::Digitize() {

    _ccdCollection = new McDarkTpcDigiCollection(GetName(), collectionName[0]); 
    _chargeCollection = new McDarkTpcDigiCollection(GetName(), collectionName[1]); 
    _pmtCollection = new McDarkTpcDigiCollection(GetName(), collectionName[2]); 

  
    G4DigiManager* DigiMan = G4DigiManager::GetDMpointer();  

    // access to TPC geometry
    G4LogicalVolume *ldet = _detector->GetLogicalVolume();
    assert(ldet);
    G4VSensitiveDetector *sdet = ldet->GetSensitiveDetector();
    assert(sdet);
    G4Tubs *tpcVolume = (G4Tubs*)ldet->GetSolid();
    double tpcZhalf = tpcVolume->GetZHalfLength();
    const G4RotationMatrix *tpcRotation=_detector->GetObjectRotation();
    G4ThreeVector tpcDriftDirection = (*tpcRotation)(G4ThreeVector(0,0,1));
    G4ThreeVector tpcXDirection = (*tpcRotation)(G4ThreeVector(1,0,0));
    G4ThreeVector tpcYDirection = (*tpcRotation)(G4ThreeVector(0,1,0));
    G4ThreeVector tpcOrigin = _detector->GetObjectTranslation()+tpcDriftDirection*tpcZhalf;


    
    // access to TPC hits
    G4int tpcID = DigiMan->GetHitsCollectionID( sdet->GetCollectionName(0) );
    McDarkTpcHitCollection *tpcHits = (McDarkTpcHitCollection*) (DigiMan->GetHitsCollection(tpcID));
    assert(tpcHits);
    
    // access to TPC gas
    McDarkTpcGas *tpcGas=(McDarkTpcGas*)ldet->GetMaterial();
    
        
    G4int nhit = tpcHits->entries();
    for (G4int i=0; i<nhit; i++) {

        McDarkTpcHit *hit = (*tpcHits)[i];

        // total number of electrons and photons after amplification
        G4double quenchedEnergy=hit->getDepositEnergy()*hit->getQuenchingFactor();
        G4double electrons=0, photons=0;
        tpcGas->energyToElectronsAndPhotons( quenchedEnergy, electrons, photons );
 
        
        // ionization point in TPC coordinates (focal plane)
        G4double dx = (tpcOrigin - hit->getDepositPosition())*tpcXDirection;
        G4double dy = (tpcOrigin - hit->getDepositPosition())*tpcYDirection;
        G4double dz = (tpcOrigin - hit->getDepositPosition())*tpcDriftDirection;
        G4ThreeVector delR( dx, dy, dz);


        //
        // detect photons with CCD and PMT
        //
        G4int chid;
        G4double weight;
        G4double photonsPerGroup=10;
        G4int nPhotonGroups=G4int(photons/photonsPerGroup);
        for (int iph=0; iph<nPhotonGroups; iph++) {
            
            G4ThreeVector tpcR = tpcGas->propagateElectrons( delR ); // local coordinates in anode plane

            /*
            G4cout << delR << "  " <<tpcR << "  " << hit->getHitID() << "   "<< hit->getTrackID()
                   << "   dE=" << hit->getDepositEnergy()/keV << "   qf="<<hit->getQuenchingFactor()
                   << "  npho="<< photons
                   << G4endl;
            */
            
            // loop over CCD's
            for (unsigned int icam=0; icam<_ccdList.size(); icam++) {
                _ccdList[icam]->processHit( tpcR, photonsPerGroup, chid, weight);
                if (chid<0) continue;
                McDarkTpcDigi* digi = new McDarkTpcDigi();
                digi->setModuleID( icam );
                digi->setChannelID( chid );
                digi->setWeight( weight );
                digi->setHitID( hit->getHitID() );
                digi->setTrackID( hit->getTrackID() );
                digi->setDigiID( _ccdCollection->entries()+1 );
                _ccdCollection->insert(digi);   
            }

            // loop over PMT's
            for (G4int ipmt=0; ipmt<0; ipmt++) {
                McDarkTpcDigi* digi = new McDarkTpcDigi();
                digi->setModuleID( ipmt );
                digi->setChannelID( -1 );
                digi->setWeight( -1 );
                digi->setHitID( hit->getHitID() );
                digi->setTrackID( hit->getTrackID() );
                digi->setDigiID( _pmtCollection->entries()+1 );
                _pmtCollection->insert(digi);   
            }

            
        }

        
        //
        // Detect electrons with charge readout
        //
        for (int iel=0; iel<electrons&&0; iel++) {
       
            for (G4int ich=0; ich<0; ich++) {
                McDarkTpcDigi* digi = new McDarkTpcDigi();
                digi->setModuleID( ich );
                digi->setChannelID( -1 );
                digi->setWeight( -1 );
                digi->setHitID( hit->getHitID() );
                digi->setTrackID( hit->getTrackID() );
                digi->setDigiID( _chargeCollection->entries()+1 );
                _chargeCollection->insert(digi);   
            }
        }

  
        
    }

    
    StoreDigiCollection(_ccdCollection);
    StoreDigiCollection(_pmtCollection);
    StoreDigiCollection(_chargeCollection);    
}



void
McDarkTpcDigitizer::addCamera(McDarkCamera *cam) { 
  _ccdList.push_back(cam); 
  
  McDarkCameraManager::addCamera(cam);
}








