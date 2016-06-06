// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// --------------------------------------------------------------
//                   15 April 1998 M.Maire
// --------------------------------------------------------------

#include "McDarkMinEkineCuts.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

McDarkMinEkineCuts::McDarkMinEkineCuts(const G4String& aName)
  : McDarkSpecialCuts(aName)
{
   if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
   SetProcessType(fUserDefined);
}

McDarkMinEkineCuts::~McDarkMinEkineCuts()
{}

McDarkMinEkineCuts::McDarkMinEkineCuts(McDarkMinEkineCuts&)
  : McDarkSpecialCuts()
{}


G4double McDarkMinEkineCuts::PostStepGetPhysicalInteractionLength(
                             const G4Track& aTrack,
			     G4double ,
			     G4ForceCondition* condition
			    )
{
  // condition is set to "Not Forced"
  *condition = NotForced;

   G4double     proposedStep = DBL_MAX;
   // get the pointer to UserLimits
   G4UserLimits* pUserLimits = aTrack.GetVolume()->GetLogicalVolume()->GetUserLimits();
   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   G4ParticleDefinition* aParticleDef = aTrack.GetDefinition();
  
   if (pUserLimits && aParticleDef->GetPDGCharge() != 0.0) {
     //min kinetic energy
     G4double temp = DBL_MAX;
     G4double    eKine     = aParticle->GetKineticEnergy();
     const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
     G4double eMin = pUserLimits->GetUserMinEkine(aTrack);

     G4double    rangeNow = DBL_MAX;

     rangeNow = G4EnergyLossTables::GetRange(aParticleDef,eKine,couple);

     if (eKine < eMin ) {
       proposedStep = 0.;
     } else {
       // charged particles only
       G4double rangeMin = G4EnergyLossTables::GetRange(aParticleDef,eMin,couple);
       temp = rangeNow - rangeMin;
       if (proposedStep > temp) proposedStep = temp;
     }
   }
   return proposedStep;
}
