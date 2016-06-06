
#include "McDarkTrackingAction.hh"
#include "McDarkRunAction.hh"
#include "McDarkCameraManager.hh"

#include "G4Track.hh"


McDarkTrackingAction::McDarkTrackingAction(McDarkRunAction* run) : _run(run)
{ }


void McDarkTrackingAction::PreUserTrackingAction(const G4Track* aTrack) {
}


void
McDarkTrackingAction::PostUserTrackingAction(const G4Track* aTrack) {

}


