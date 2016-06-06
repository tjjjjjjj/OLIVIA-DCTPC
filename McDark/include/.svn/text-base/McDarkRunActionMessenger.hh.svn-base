//
// $Id: McDarkRunActionMessenger.hh,v 1.4 2009/12/30 22:34:49 ddujmic Exp $
// GEANT4 tag $Name:  $
//
// --------------------------------------------------------------


#ifndef McDarkRunActionMessenger_h
#define McDarkRunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class McDarkRunAction;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class McDarkRunActionMessenger: public G4UImessenger
{
  public:
    McDarkRunActionMessenger(McDarkRunAction*);
   ~McDarkRunActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    McDarkRunAction*   McDarkRun;


    G4UIcmdWithABool*     PlotEventCmd;
    G4UIcmdWithAnInteger*     PrintTracksCmd;
    G4UIcmdWithABool*     PrintHitsCmd;
};

#endif

