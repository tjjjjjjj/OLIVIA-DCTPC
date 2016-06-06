// $Id: McDarkSteppingAction.cc,v 1.11 2009/12/30 22:41:07 ddujmic Exp $
// GEANT4 tag $Name:  $
// 

#include "McDarkSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "McDarkRunAction.hh"
#include "G4RunManager.hh"
#include "globals.hh"



McDarkSteppingAction::McDarkSteppingAction(McDarkRunAction *run) : _run(run)
{}

void McDarkSteppingAction::UserSteppingAction(const G4Step* aStep) {
    return;
}
