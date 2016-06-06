
#include "McDarkCameraManagerMessenger.hh"
#include "McDarkCameraManager.hh"
#include "McDarkUICameraCmd.hh"

#include <sstream>

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "McDarkUIcmdWithAnIntegerAndDouble.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"

//#include "McDarkSpergelDistribution.hh"

McDarkCameraManagerMessenger::McDarkCameraManagerMessenger(McDarkCameraManager* camMan)
: McDarkCameraMan(camMan) {
    
}


McDarkCameraManagerMessenger::~McDarkCameraManagerMessenger() {

}


void McDarkCameraManagerMessenger::SetNewValue(G4UIcommand * command,G4String newValue) {

}
