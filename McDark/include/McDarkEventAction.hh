// $Id: McDarkEventAction.hh,v 1.9 2010/01/04 21:49:30 ddujmic Exp $
// GEANT4 tag $Name:  $
//
 
#ifndef McDarkEventAction_h
#define McDarkEventAction_h 1

#include "G4UserEventAction.hh"
#include "McDarkCameraManager.hh"
#include <iostream>

#include "TH2F.h"

class G4Event;
class McDarkRunAction;
//class McDarkCameraManager;
class G4TrajectoryContainer;

class McDarkEventAction : public G4UserEventAction {
    
  public:
    McDarkEventAction(McDarkRunAction *run);
    ~McDarkEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

    McDarkRunAction* getRun() { return _run; }

protected:
    
    void printTrackTable(const G4Event* evt);
    
    void exportTracks(const G4Event* evt);
    void exportCcdInfo(const G4Event* evt);
    
    G4bool isTriggered(const G4Event* evt);
    
private:
    McDarkRunAction* _run;
    McDarkCameraManager* camMan;
};


#endif

    
