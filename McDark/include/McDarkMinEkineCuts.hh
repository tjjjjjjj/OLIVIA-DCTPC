// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// ------------------------------------------------------------
//                  14 Aug. 1998  H.Kurashige
// ------------------------------------------------------------

#ifndef DMXMinEkineCuts_h
#define DMXMinEkineCuts_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "McDarkSpecialCuts.hh"


class McDarkMinEkineCuts : public McDarkSpecialCuts
{
  public:     

     McDarkMinEkineCuts(const G4String& processName ="McDarkMinEkineCuts" );

     virtual ~McDarkMinEkineCuts();

     // PostStep GPIL
     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );
            
			    
  private:
  
  // hide assignment operator as private 
      McDarkMinEkineCuts(McDarkMinEkineCuts&);
      McDarkMinEkineCuts& operator=(const McDarkMinEkineCuts& right);

};

#endif

