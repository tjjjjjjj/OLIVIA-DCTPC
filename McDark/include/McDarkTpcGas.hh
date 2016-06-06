#ifndef McDarkTpcGas_h
#define McDarkTpcGas_h 1

#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4ParticleDefinition;


class McDarkTpcGas : public G4Material {
    
public:

    McDarkTpcGas(const G4String& name,                              //its name
                 G4double  density,                           //density
                 G4int     nComponents,                       //nbOfComponents
                 G4State   state    = kStateUndefined,        //solid,gas
                 G4double  temp     = STP_Temperature,        //temperature
                 G4double  pressure = STP_Pressure);          //pressure
    
    void setWorkFunction(G4double energy) { _workFunction=energy; }
    void setAvalancheGain(G4double gain) { _avalancheGain=gain; }
    void setPhotonsPerElectron(G4double ppe) { _photonsPerElectron=ppe; }
    void setTDiffusion(G4double coef) { _tDiffCoef=coef; }
    void setLDiffusion(G4double coef) { _lDiffCoef=coef; }
    void setTSpread(G4double spread) { _tSpreadSq=spread*spread; }
    void setLSpread(G4double spread) { _lSpreadSq=spread*spread; }
    void setDriftV(G4double v) { _driftV=v; }

    void energyToElectronsAndPhotons(G4double quenchedEnergy, G4double &electrons, G4double &photons);

    G4ThreeVector propagateElectrons(G4ThreeVector delrad);
    
protected:
    
    G4double energyToElectrons(G4double quenchedEnergy);
    G4double electronsToPhotons(G4double electrons);
    
    G4double smearX(G4double dx, G4double dz);
    G4double smearY(G4double dy, G4double dz) { return smearX(dy, dz); }
    G4double smearZ(G4double  dz);

    
private:
  
    G4double _workFunction; // ionization work function
    G4double _avalancheGain; // avalanche amplification gain
    G4double _photonsPerElectron; // number of scintillation photons per ionization electron
    G4double _tDiffCoef; // transverse diffusion
    G4double _lDiffCoef; // longitudinal diffusion coefficient 
    G4double _tSpreadSq; // spread due to mesh pitch, avalanche
    G4double _lSpreadSq; // spread due to mesh pitch, avalanche
    G4double _driftV; // drift velocity in sensitive volume
};




#endif








