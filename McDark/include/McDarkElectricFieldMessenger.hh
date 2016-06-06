//
// Created by: T. Sahin (tcsahin@MIT.EDU), D. Dujmic (ddujmic@MIT.EDU)
// Date:       Set 15, 2007
// Copyright:  MIT 2007
//
//
// $Id: McDarkElectricFieldMessenger.hh,v 1.1 2007/10/17 11:51:17 ddujmic Exp $
// GEANT4 tag $Name:  $
//
// 


#ifndef McDarkElectricFieldMessenger_h
#define McDarkElectricFieldMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class McDarkElectricFieldSetup;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;


class McDarkElectricFieldMessenger: public G4UImessenger
{
  public:
    McDarkElectricFieldMessenger(McDarkElectricFieldSetup* );
   ~McDarkElectricFieldMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    void SetNewValue(G4UIcommand*, G4int);
    
  private:

    McDarkElectricFieldSetup*  fEFieldSetup;
    
    G4UIdirectory*             mcDarkdetDir;

    G4UIcmdWithAnInteger*      stepperCmd;
    G4UIcmdWithADoubleAndUnit* elecFieldCmd;
    G4UIcmdWithADoubleAndUnit* minStepCmd;
    G4UIcmdWithoutParameter*   updateCmd;
    G4UIcmdWithAString*        dielectricCmd;
    G4UIcmdWithADoubleAndUnit* cathodeVoltageCmd;
    G4UIcmdWithADoubleAndUnit* anodeVoltageCmd;

};

#endif

