#include "McDarkTpcSD.hh"
#include "G4VPhysicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Ions.hh"
#include "G4ios.hh"

#include "G4hNuclearStoppingModel.hh"
#include "G4hParametrisedLossModel.hh"


McDarkTpcSD::McDarkTpcSD(G4String name) : G4VSensitiveDetector(name) {
  collectionName.insert( name+"_hitCollection" );
}


McDarkTpcSD::~McDarkTpcSD(){ }


void
McDarkTpcSD::Initialize(G4HCofThisEvent*) {
  _hitCollection = new McDarkTpcHitCollection( SensitiveDetectorName, collectionName[0]);
}


G4bool
McDarkTpcSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    
  G4double edep = aStep->GetTotalEnergyDeposit();
  if ( edep<50*eV ) return false;

  
  // track info
  G4Track *track=aStep->GetTrack();
  G4ParticleDefinition* particle = track->GetDefinition();
  if ( fabs(particle->GetPDGCharge())<1e-6 ) return false;
   
  
  // hit info
  McDarkTpcHit* newHit = new McDarkTpcHit();
  newHit->setDepositEnergy(edep);
  G4ThreeVector pos = 0.5 * ( aStep->GetPreStepPoint()->GetPosition() +
                              aStep->GetPostStepPoint()->GetPosition() );
  newHit->setDepositPosition(pos);
  double time = 0.5 * ( aStep->GetPreStepPoint()->GetGlobalTime() +
                        aStep->GetPostStepPoint()->GetGlobalTime() );
  newHit->setDepositTime(time);
  
  // track index
  newHit->setTrackID( track->GetTrackID() );

  // hit index
  newHit->setHitID( _hitCollection->entries() );

  // track PDG
  newHit->setPDGIndex( particle->GetPDGEncoding() );

  // quenching factor
  newHit->setQuenchingFactor( 0.0 );

  if (particle->GetPDGCharge()!=0.0) switch (particle->GetPDGEncoding()) {

  // for charged leptons assume all energy is 'electron equivalent'
  // ... note: if a nucleus creates ionization electrons, both will create hits. 
  // But most of ioniation electrons fall below the energy threshold so the hits
  // are created mostly from primary ionization particle and overcounting of
  // deposited energy is negligible  

  case 11:
  case 13:
    newHit->setQuenchingFactor( 1.0 );
    break;

  // for nuclei find quenching factor relative to electrons 
  default:
    G4TouchableHistory* hist = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    assert(hist);
    const G4VPhysicalVolume* physVol = hist->GetVolume();
    assert(physVol);
    G4LogicalVolume *logVol = physVol->GetLogicalVolume();
    assert(logVol);
    G4Material *volMaterial=logVol->GetMaterial();
    assert(volMaterial);
    G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();
    static G4hParametrisedLossModel* theElectronicStoppingModel= new G4hParametrisedLossModel("SRIM2000p");
    G4double Se = theElectronicStoppingModel->TheValue(particle, volMaterial, energy);
    static G4hNuclearStoppingModel* theNuclearStoppingModel = new G4hNuclearStoppingModel("ICRU_R49") ;
    G4double Sn = theNuclearStoppingModel->TheValue(particle, volMaterial, energy);
    newHit->setQuenchingFactor( (Sn*0.3 + Se)/(Sn + Se) );
    break;
  }


  
  _hitCollection->insert(newHit);
  
  return true;
}


void
McDarkTpcSD::EndOfEvent(G4HCofThisEvent* HCE) {

  G4String HCname = collectionName[0];
  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  G4int HCID = SDman->GetCollectionID(HCname);
  HCE->AddHitsCollection(HCID, _hitCollection);
  //G4int nHits = _hitCollection->entries();
  //if (nHits)  G4cout << GetName() << ": Collection " << HCname << " with ID " << HCID << " has " <<  nHits << " hits" << G4endl;
}


void
McDarkTpcSD::clear() {} 


void
McDarkTpcSD::DrawAll() {} 


void
McDarkTpcSD::PrintAll() {} 


