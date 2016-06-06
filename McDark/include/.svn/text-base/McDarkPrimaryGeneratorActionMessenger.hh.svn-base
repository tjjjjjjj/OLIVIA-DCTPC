//
// $Id: McDarkPrimaryGeneratorActionMessenger.hh,v 1.4 2009/11/28 23:12:31 ddujmic Exp $
// GEANT4 tag $Name:  $
//
// --------------------------------------------------------------


#ifndef McDarkPrimaryGeneratorActionMessenger_h
#define McDarkPrimaryGeneratorActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class McDarkPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIdirectory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class McDarkPrimaryGeneratorActionMessenger: public G4UImessenger
{
  public:
    McDarkPrimaryGeneratorActionMessenger(McDarkPrimaryGeneratorAction*);
   ~McDarkPrimaryGeneratorActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    McDarkPrimaryGeneratorAction*   McDarkPrimaryGenerator;

    G4UIcmdWithADouble*   wimpMassCmd;
    G4UIcmdWithADouble*   wimpDensityCmd;
    G4UIcmdWithADouble*   mwVelocityCmd;
    G4UIcmdWithADouble*   earthVelocityCmd;
    G4UIcmdWithADouble*   earthDVelocityCmd;
    G4UIcmdWithADouble*   yearFractionCmd;

    G4UIcmdWithAString*   distributionCmd;

    G4UIcmdWithADoubleAndUnit*	minEnergyCmd;
    G4UIcmdWithADoubleAndUnit*  maxEnergyCmd;
    G4UIcmdWithADoubleAndUnit*	minXposCmd;
    G4UIcmdWithADoubleAndUnit*  maxXposCmd;
    G4UIcmdWithADoubleAndUnit*  minYposCmd;
    G4UIcmdWithADoubleAndUnit*  maxYposCmd;
    G4UIcmdWithADoubleAndUnit*  minZposCmd;
    G4UIcmdWithADoubleAndUnit*	maxZposCmd;


    G4UIdirectory*        generatorDirectory;
    G4UIdirectory*	  spergelDirectory;
    G4UIdirectory*	  isotropicDirectory;

};

#endif

