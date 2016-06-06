//
// $Id: McDarkIonMultipleScattering.hh,v 1.1 2008/03/11 20:44:22 ddujmic Exp $
// GEANT4 tag $Name:  $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef McDarkIonMultipleScattering_h
#define McDarkIonMultipleScattering_h 1

#include "G4VMultipleScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class McDarkUrbanMscModel;

class McDarkIonMultipleScattering : public G4VMultipleScattering

{
public:    // with description

  McDarkIonMultipleScattering(const G4String& processName="msc");

  virtual ~McDarkIonMultipleScattering();

  // returns true for charged particles, false otherwise
  G4bool IsApplicable (const G4ParticleDefinition& p);

  // Print few lines of informations about the process: validity range,
  void PrintInfo();

  // set boolean flag steppingAlgorithm
  // ( true/false : standard or 7.1 style process)
  void MscStepLimitation(G4bool algorithm, G4double factor = -1.);

  // geom. step length distribution should be sampled or not
  void Setsamplez(G4bool value) { samplez = value;};

  // to reduce the energy/step dependence
  void Setdtrl(G4double value) { dtrl = value;};

  // 'soften' step limitation above lambdalimit
  void SetLambdalimit(G4double value) { lambdalimit = value;};

  // Steplimit = facrange*max(range,lambda)
  void SetFacrange(G4double val) { facrange=val;};

  // connected with step size reduction due to geometry
  void SetFacgeom(G4double val) { facgeom=val;};

  // set msc parameter skin
  // if skin <= 0 --> no single scattering at boundary
  void SetSkin(G4double val) { skin=val;};

protected:

  // This function initialise models
  void InitialiseProcess(const G4ParticleDefinition*);

private:        // data members

  McDarkUrbanMscModel* mscUrban;

  G4double lowKineticEnergy;
  G4double highKineticEnergy;
  G4int    totBins;

  G4double lambdalimit;
  G4double facrange;
  G4double facgeom;
  G4double skin; 
  G4double dtrl;

  G4bool   steppingAlgorithm;
  G4bool   samplez;
  G4bool   isInitialized;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
