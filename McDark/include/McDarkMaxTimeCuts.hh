// Stolen from:
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// MaxTimeCuts header
// --------------------------------------------------------------

#ifndef McDarkMaxTimeCuts_h
#define McDarkMaxTimeCuts_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "McDarkSpecialCuts.hh"


class McDarkMaxTimeCuts : public McDarkSpecialCuts
{
  public:     

     McDarkMaxTimeCuts(const G4String& processName ="McDarkMaxTimeCuts" );

     virtual ~McDarkMaxTimeCuts();

     // PostStep GPIL
     virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );
            
			    
  private:
  
  // hide assignment operator as private 
      McDarkMaxTimeCuts(McDarkMaxTimeCuts&);
      McDarkMaxTimeCuts& operator=(const McDarkMaxTimeCuts& right);

};

#endif

