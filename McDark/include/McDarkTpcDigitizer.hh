#ifndef McDarkTpcDigitizer_h
#define McDarkTpcDigitizer_h 1

#include "G4VDigitizerModule.hh"
#include "G4VUserDetectorConstruction.hh"
#include "McDarkTpcDigi.hh"
#include "globals.hh"

class McDarkCamera;

#include <vector>
using std::vector;

class G4VSensitiveDetector;
class G4PVPlacement;


class McDarkTpcDigitizer : public G4VDigitizerModule {
    
public:
  
    McDarkTpcDigitizer(G4String name, const G4PVPlacement *det);
    
    ~McDarkTpcDigitizer();
    
    void Digitize();

    void addCamera(McDarkCamera *cam); 

    McDarkCamera* ccd(unsigned int icam) { return _ccdList[icam]; }

    unsigned int nccd() { return _ccdList.size(); }

    McDarkTpcDigi* ccdDigi(int i) { return (*_ccdCollection)[i]; }
    int nDigiCCD() { return _ccdCollection->entries(); }

protected:

    
private:
  
    McDarkTpcDigiCollection*  _ccdCollection; // pointer to digi container

    McDarkTpcDigiCollection*  _chargeCollection; // pointer to digi container

    McDarkTpcDigiCollection*  _pmtCollection; // pointer to digi container

    const G4PVPlacement *_detector; // pointer to TPC

    vector<McDarkCamera*> _ccdList;

    
};



#endif








