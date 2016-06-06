#include "McDarkWormDigitizer.hh"
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
#include "G4Box.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleTable.hh"

McDarkWormDigitizer::McDarkWormDigitizer(G4String name,  McDarkCamera *cam)
    : G4VDigitizerModule(name),
      _cam(cam) {

    McDarkCameraManager::addCamera(cam);

    G4String colName = name+"_ccdCollection";
    collectionName.push_back(colName); 
}


McDarkWormDigitizer::~McDarkWormDigitizer() {}


void
McDarkWormDigitizer::Digitize() {
  
    _ccdCollection = new McDarkTpcDigiCollection(GetName(), collectionName[0]); 

  
    G4DigiManager* DigiMan = G4DigiManager::GetDMpointer();  

    // access to TPC geometry
    G4LogicalVolume *ldet = _cam->getActiveVolume()->GetLogicalVolume();
    assert(ldet);
    G4VSensitiveDetector *sdet = ldet->GetSensitiveDetector();
    assert(sdet);
    G4Box *ccdVolume = (G4Box*)ldet->GetSolid();
    double ccdXhalf = ccdVolume->GetXHalfLength();
    double ccdYhalf = ccdVolume->GetYHalfLength();
    double ccdZhalf = ccdVolume->GetZHalfLength();
    const G4RotationMatrix *ccdRotation=_cam->getActiveVolume()->GetObjectRotation();
    G4ThreeVector ccdDriftDirection = (*ccdRotation)(G4ThreeVector(0,0,1));
    G4ThreeVector ccdXDirection = (*ccdRotation)(G4ThreeVector(1,0,0));
    G4ThreeVector ccdYDirection = (*ccdRotation)(G4ThreeVector(0,1,0));
    G4ThreeVector ccdOrigin = _cam->getActiveVolume()->GetObjectTranslation()+ccdDriftDirection*ccdZhalf;

    
    // access to CCD worms
    G4int ccdID = DigiMan->GetHitsCollectionID( sdet->GetCollectionName(0) );
    McDarkTpcHitCollection *ccdHits = (McDarkTpcHitCollection*) (DigiMan->GetHitsCollection(ccdID));
    assert(ccdHits);
    
        
    G4int nhit = ccdHits->entries();
    for (G4int i=0; i<nhit; i++) {

        McDarkTpcHit *hit = (*ccdHits)[i];

        // total number of electrons and photons after amplification
        G4double quenchedEnergy=hit->getDepositEnergy()*hit->getQuenchingFactor();
        
        // ionization point in TPC coordinates (focal plane)
        G4double dx = (ccdOrigin - hit->getDepositPosition())*ccdXDirection/(2*ccdXhalf)+0.5;
        G4double dy = (ccdOrigin - hit->getDepositPosition())*ccdYDirection/(2*ccdYhalf)+0.5;
        G4double dz = 0;
        G4ThreeVector ccdR( dx, dy, dz);


        //
        // detect photons with CCD and PMT
        //
        G4int chid;
        G4double weight;            
	_cam->processWorm( ccdR, quenchedEnergy, chid, weight);
	if (chid<0) continue;
	McDarkTpcDigi* digi = new McDarkTpcDigi();
	digi->setModuleID( 0 );
	digi->setChannelID( chid );
	digi->setWeight( weight );
	digi->setHitID( hit->getHitID() );
	digi->setTrackID( hit->getTrackID() );
	digi->setDigiID( _ccdCollection->entries()+1 );
	_ccdCollection->insert(digi);        
    }    

    
    StoreDigiCollection(_ccdCollection);    
}









