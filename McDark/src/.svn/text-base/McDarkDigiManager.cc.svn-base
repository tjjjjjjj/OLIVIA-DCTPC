#include "McDarkDigiManager.hh"


#include "assert.h"
#include<vector>
using std::vector;

vector<G4String> McDarkDigiManager::listOfDigiModules;


G4DigiManager* 
McDarkDigiManager::getDMpointer() { 
  return G4DigiManager::GetDMpointer(); 
}


void 
McDarkDigiManager::Digitize() {
  
  for (unsigned int i=0; i<listOfDigiModules.size(); i++){
    getDMpointer()->Digitize( listOfDigiModules[i] );
  }

}


void 
McDarkDigiManager::addNewModule(G4VDigitizerModule* DM) {
  
  getDMpointer()->AddNewModule( DM );
  listOfDigiModules.push_back( DM->GetName() );  

}
