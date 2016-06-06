#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "McDarkDetectorConstruction.hh"
#include "McDarkNeutronRadiography.hh"
#include "McDarkPhysicsList.hh"
#include "McDarkPrimaryGeneratorAction.hh"

#include "McDarkRunAction.hh"
#include "McDarkEventAction.hh"
#include "McDarkSteppingAction.hh"
#include "McDarkTrackingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE_ROOT
#include "G4UIRoot.hh"
#endif

int main(int argc, char** argv) {
    
  // Construct the default run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  //
  G4VUserDetectorConstruction* detector = new McDarkDetectorConstruction;
  //G4VUserDetectorConstruction* detector = new McDarkNeutronRadiography;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new McDarkPhysicsList;
  runManager->SetUserInitialization(physics);
  
  // set mandatory user action class
  //
  G4VUserPrimaryGeneratorAction* gen_action = new McDarkPrimaryGeneratorAction;
  runManager->SetUserAction(gen_action);

  G4UserRunAction* run_action = new McDarkRunAction;
  runManager->SetUserAction(run_action);

  G4UserEventAction* event_action = new McDarkEventAction((McDarkRunAction*)run_action);
  runManager->SetUserAction(event_action);

  G4UserSteppingAction* stepping_action = new McDarkSteppingAction((McDarkRunAction*)run_action);
  runManager->SetUserAction(stepping_action);

  G4UserTrackingAction* tracking_action = new McDarkTrackingAction((McDarkRunAction*)run_action);
  runManager->SetUserAction(tracking_action);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // Initialize G4 kernel
  //
  runManager->Initialize();

  // Get the pointer to the UI manager and set verbosities
  //
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if(argc==1)  // Define (G)UI terminal for interactive mode
  { 
    // G4UIterminal is a (dumb) terminal
    //
    G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif
#ifdef G4UI_USE_ROOT
  // session = new G4UIRoot(argc, argv);
#endif
    UI->ApplyCommand("/control/execute vis.mac");    
    session->SessionStart();
    delete session;
  }
  else   // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

  // Start a run
  //
  // G4int numberOfEvent = 3;
  // runManager->BeamOn(numberOfEvent);

  // Job termination
  //
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  //
  
#ifdef G4VIS_USE
  delete visManager;
#endif

  
  delete runManager;

  return 0;
}




