
#ifndef McDarkCameraManagerMessenger_h
#define McDarkCameraManagerMessenger_h 1

// --------------------------------------------------------------
#include "globals.hh"
#include "G4UImessenger.hh"

class McDarkCameraManager;
class McDarkUICameraCmd;
class McDarkUIcmdWithAnIntegerAndDouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIdirectory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class McDarkCameraManagerMessenger: public G4UImessenger
{
  public:
    McDarkCameraManagerMessenger(McDarkCameraManager*);
   ~McDarkCameraManagerMessenger();

    void SetNewValue(G4UIcommand *command, G4String newValue);
    
  private:
    McDarkCameraManager*   McDarkCameraMan;

    G4UIcmdWithAnInteger*	addCameraCmd;
    McDarkUICameraCmd*		setupCameraCmd;
    McDarkUIcmdWithAnIntegerAndDouble* setNoiseCmd;
    G4UIcmdWithABool*		toggleNoiseCmd;

    G4UIdirectory*	  cameraDirectory;

};

#endif

