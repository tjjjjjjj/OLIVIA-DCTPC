// $Id: McDarkSteppingAction.hh,v 1.8 2009/12/30 22:41:07 ddujmic Exp $
// GEANT4 tag $Name:  $
// 

#ifndef McDarkSteppingAction_h
#define McDarkSteppingAction_h 

#include "G4UserSteppingAction.hh"
#include "globals.hh"


class McDarkCameraManager;
class McDarkRunAction;
class G4ParticleDefinition;
class G4Step;

class McDarkSteppingAction : public G4UserSteppingAction {
    
  public:
    McDarkSteppingAction(McDarkRunAction *run);
   ~McDarkSteppingAction(){};

    void UserSteppingAction(const G4Step*);

    McDarkRunAction* getRun() { return _run; }
    
private:
    McDarkRunAction* _run;
};


#endif
