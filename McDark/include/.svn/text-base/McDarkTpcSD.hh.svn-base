#ifndef McDarkTpcSD_h
#define McDarkTpcSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "McDarkTpcHit.hh"

class G4Step;
class G4HCofThisEvent;



class McDarkTpcSD : public G4VSensitiveDetector {
    
  public:
  
    McDarkTpcSD(G4String);
    ~McDarkTpcSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);
    void clear();
    void DrawAll();
    void PrintAll();
  
private:
  
    McDarkTpcHitCollection*  _hitCollection;

    double quenchingFactor;
    double transverseDiffusion;
    double longitudinalDiffusion;
    
};

#endif

