#ifndef McDarkDigiManager_H
#define McDarkDigiManager_H 1

#include "globals.hh"
#include "G4String.hh"
#include "G4DigiManager.hh"

#include <vector>
using std::vector;

class McDarkDigiManager {

public:


  static void Digitize();

  static G4DigiManager* getDMpointer();

  static vector<G4String> listOfDigiModules;

  static void addNewModule(G4VDigitizerModule* DM);

};
    
#endif
