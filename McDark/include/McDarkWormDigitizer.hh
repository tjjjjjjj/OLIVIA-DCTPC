#ifndef McDarkWormDigitizer_h
#define McDarkWormDigitizer_h 1

#include "G4VDigitizerModule.hh"
#include "G4VUserDetectorConstruction.hh"
#include "McDarkTpcDigi.hh"
#include "globals.hh"

class McDarkCamera;

#include <vector>
using std::vector;

class G4VSensitiveDetector;
class G4PVPlacement;
class McDarkCamera;

class McDarkWormDigitizer : public G4VDigitizerModule {
    
public:
  
    McDarkWormDigitizer(G4String name,  McDarkCamera *cam);
    
    ~McDarkWormDigitizer();
    
    void Digitize();

    McDarkCamera* ccd() { return _cam; }

    McDarkTpcDigi* ccdDigi(int i) { return (*_ccdCollection)[i]; }
    int nDigiCCD() { return _ccdCollection->entries(); }

protected:


    
private:
  
    McDarkTpcDigiCollection*  _ccdCollection; // pointer to digi container

    McDarkCamera* _cam;
};



#endif








