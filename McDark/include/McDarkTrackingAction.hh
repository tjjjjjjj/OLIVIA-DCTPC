#ifndef McDarkTrackingAction_h
#define McDarkTrackingAction_h

#include "G4UserTrackingAction.hh"

class McDarkCameraManager;
class McDarkRunAction;
class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class McDarkTrackingAction : public G4UserTrackingAction {

  public:  
    McDarkTrackingAction(McDarkRunAction* );
    
   ~McDarkTrackingAction() {};
   
    void  PreUserTrackingAction(const G4Track*);
    
    void PostUserTrackingAction(const G4Track*);
         
    McDarkRunAction* getRun() { return _run; }
    
  private:
    McDarkRunAction*              _run;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
