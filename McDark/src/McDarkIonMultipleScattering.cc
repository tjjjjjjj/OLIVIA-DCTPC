//
// $Id: McDarkIonMultipleScattering.cc,v 1.3 2010/12/29 13:14:53 ddujmic Exp $
// GEANT4 tag $Name:  $


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "McDarkIonMultipleScattering.hh"
#include "McDarkUrbanMscModel.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

McDarkIonMultipleScattering::McDarkIonMultipleScattering(const G4String& processName)
  : G4VMultipleScattering(processName)
{
  lowKineticEnergy  = 100*eV;
  highKineticEnergy = 100.*MeV;
  totBins           = 12000;

  facrange          = 0.002;
  dtrl              = 0.005;
  lambdalimit       = 0.001*mm;
  facgeom           = 0.005;
  // there is no single scattering for this skin <= 0  
  // to have single scattering at boundary 
  //  skin should be > 0 ! 
  skin              = 0.;
  
  steppingAlgorithm = false;
  samplez           = false ; 
  isInitialized     = false;  

  SetBinning(totBins);
  SetMinKinEnergy(lowKineticEnergy);
  SetMaxKinEnergy(highKineticEnergy);

  SetLateralDisplasmentFlag(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

McDarkIonMultipleScattering::~McDarkIonMultipleScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool McDarkIonMultipleScattering::IsApplicable (const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void McDarkIonMultipleScattering::MscStepLimitation(G4bool algorithm, G4double factor) 
{
  steppingAlgorithm = algorithm;
  if (factor > 0.) SetFacrange(factor);
  else { if (algorithm) SetFacrange(0.02); else SetFacrange(0.2);}

  if(verboseLevel > 1)  
    G4cout << "Stepping algorithm is set to " << steppingAlgorithm 
	   << " with facrange = " << facrange << G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void McDarkIonMultipleScattering::InitialiseProcess(const G4ParticleDefinition* p)
{
  if(isInitialized) {
    //mscUrban->SetMscStepLimitation(steppingAlgorithm, facrange);
    return;
  }

//   if (p->GetParticleType() == "nucleus") {
//     SetLateralDisplasmentFlag(false);
//     SetBuildLambdaTable(false);
//   } else {
//     SetBuildLambdaTable(true);
//   }
  /*SetLateralDisplasmentFlag(true);//dd
  SetBuildLambdaTable(true);//dd
  mscUrban = new McDarkUrbanMscModel(facrange,dtrl,lambdalimit,
				     facgeom,skin,
				     samplez,steppingAlgorithm);
  mscUrban->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
  mscUrban->SetLowEnergyLimit(lowKineticEnergy);
  mscUrban->SetHighEnergyLimit(highKineticEnergy);
  AddEmModel(1,mscUrban); 
  isInitialized = true; */

  mscUrban = new McDarkUrbanMscModel();
  AddEmModel(1,mscUrban);
  isInitialized = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void McDarkIonMultipleScattering::PrintInfo()
{
  G4cout << "      Boundary/stepping algorithm is active with facrange= "
	 << facrange
	 << "  Step limitation " << steppingAlgorithm
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

