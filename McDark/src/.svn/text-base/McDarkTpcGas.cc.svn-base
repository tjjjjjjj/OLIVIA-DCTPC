#include "McDarkTpcGas.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

McDarkTpcGas::McDarkTpcGas(const G4String& name,                           
                           G4double  density,                          
                           G4int     nComponents,                      
                           G4State   state,       
                           G4double  temp,        
                           G4double  pressure) :         
    G4Material(name, density, nComponents, state, temp, pressure),
    _workFunction(34.3*eV),
    _avalancheGain(1e5),
    _photonsPerElectron(0.0),
    _tDiffCoef(0.0),
    _lDiffCoef(0.0),
    _tSpreadSq(0.0),
    _lSpreadSq(0.0),
    _driftV(1e7*cm/s)
{}
 

void
McDarkTpcGas::energyToElectronsAndPhotons(G4double energy, G4double &electrons, G4double &photons) {
    // Converts energy to number of electrons on the anode and
    // number of scintillation photons
    // 
    
    G4double meanElectrons = energyToElectrons(energy);
    G4double meanPhotons = electronsToPhotons(meanElectrons);
    
    electrons = CLHEP::RandPoissonQ::shoot(CLHEP::HepRandom::getTheEngine(), meanElectrons); 
    photons = CLHEP::RandPoissonQ::shoot(CLHEP::HepRandom::getTheEngine(), meanPhotons);
}


G4double
McDarkTpcGas::energyToElectrons(G4double quenchedEnergy) {
    return quenchedEnergy/_workFunction*_avalancheGain;
}


G4double
McDarkTpcGas::electronsToPhotons(G4double electrons) {
    return electrons*_photonsPerElectron;
}


G4double
McDarkTpcGas::smearX(G4double dx, G4double dz) {
    // assumes mean spread is dx at anode; make random spread around this value
    G4double totSpread = sqrt( _tSpreadSq + _tDiffCoef*dz );
    return CLHEP::RandGaussQ::shoot(CLHEP::HepRandom::getTheEngine(), dx, totSpread);
}


G4double
McDarkTpcGas::smearZ(G4double dz) {
    // assumes mean dz=0 at anode; make random spread around this value
    G4double totSpread = sqrt( _lSpreadSq + _lDiffCoef*dz );
    return CLHEP::RandGaussQ::shoot(CLHEP::HepRandom::getTheEngine(), 0.0, totSpread);   
}

G4ThreeVector
McDarkTpcGas::propagateElectrons(G4ThreeVector dR) {
    return G4ThreeVector( smearX(dR.x(),dR.z()),   smearY(dR.y(),dR.z()),  smearZ(dR.z()) );
}
