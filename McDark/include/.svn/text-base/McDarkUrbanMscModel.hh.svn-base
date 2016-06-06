// $Id: McDarkUrbanMscModel.hh,v 1.3 2010/12/29 13:15:06 ddujmic Exp $
// GEANT4 tag $Name:  $
//
// -------------------------------------------------------------------
//
//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// H.W.Lewis Phys Rev 78 (1950) 526 and L.Urban model

// -------------------------------------------------------------------
//

#ifndef G4UrbanMscModel93_h
#define McDarkUrbanMscModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4VMscModel.hh"
#include "G4PhysicsTable.hh"
#include "G4MscStepLimitType.hh"

class G4ParticleChangeForMSC;
class G4SafetyHelper;
class G4LossTableManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class McDarkUrbanMscModel : public G4VMscModel
{

public:

  McDarkUrbanMscModel(const G4String& nam = "UrbanMsc93");

  virtual ~McDarkUrbanMscModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition* particle,
				      G4double KineticEnergy,
				      G4double AtomicNumber,
				      G4double AtomicWeight=0., 
				      G4double cut =0.,
				      G4double emax=DBL_MAX);

  void SampleScattering(const G4DynamicParticle*,
			G4double safety);

  G4double ComputeTruePathLengthLimit(const G4Track& track,
				      G4PhysicsTable* theLambdaTable,
				      G4double currentMinimalStep);

  G4double ComputeGeomPathLength(G4double truePathLength);

  G4double ComputeTrueStepLength(G4double geomStepLength);

  G4double ComputeTheta0(G4double truePathLength,
                         G4double KineticEnergy);

private:

  G4double SimpleScattering(G4double xmeanth, G4double x2meanth);

  G4double SampleCosineTheta(G4double trueStepLength, G4double KineticEnergy);

  G4double SampleDisplacement();

  G4double LatCorrelation();

  inline G4double GetLambda(G4double kinEnergy);

  inline void SetParticle(const G4ParticleDefinition*);

  inline void UpdateCache();

  //  hide assignment operator
  McDarkUrbanMscModel & operator=(const  McDarkUrbanMscModel &right);
  McDarkUrbanMscModel(const  McDarkUrbanMscModel&);

  const G4ParticleDefinition* particle;
  G4ParticleChangeForMSC*     fParticleChange;

  G4PhysicsTable*             theLambdaTable;
  const G4MaterialCutsCouple* couple;
  G4LossTableManager*         theManager;

  G4double mass;
  G4double charge,ChargeSquare;
  G4double masslimite,lambdalimit,fr;

  G4double taubig;
  G4double tausmall;
  G4double taulim;
  G4double currentTau;
  G4double tlimit;
  G4double tlimitmin;
  G4double tlimitminfix;
  G4double tgeom;

  G4double geombig;
  G4double geommin;
  G4double geomlimit;
  G4double skindepth;
  G4double smallstep;

  G4double presafety;

  G4double lambda0;
  G4double lambdaeff;
  G4double tPathLength;
  G4double zPathLength;
  G4double par1,par2,par3;

  G4double stepmin;

  G4double currentKinEnergy;
  G4double currentRange; 
  G4double rangeinit;
  G4double currentRadLength;

  G4double theta0max,rellossmax;
  G4double third;

  G4int    currentMaterialIndex;

  G4double y;
  G4double Zold;
  G4double Zeff,Z2,Z23,lnZ;
  G4double coeffth1,coeffth2;
  G4double coeffc1,coeffc2;
  G4double scr1ini,scr2ini,scr1,scr2;

  G4bool   isInitialized;
  G4bool   inside;
  G4bool   insideskin;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4double McDarkUrbanMscModel::GetLambda(G4double e)
{
  G4double x;
  if(theLambdaTable) {
    G4bool b;
    x = ((*theLambdaTable)[currentMaterialIndex])->GetValue(e, b);
  } else {
    x = CrossSection(couple,particle,e);
  }
  if(x > DBL_MIN) x = 1./x;
  else            x = DBL_MAX;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void McDarkUrbanMscModel::SetParticle(const G4ParticleDefinition* p)
{
  if (p != particle) {
    particle = p;
    mass = p->GetPDGMass();
    charge = p->GetPDGCharge()/eplus;
    ChargeSquare = charge*charge;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void McDarkUrbanMscModel::UpdateCache()                                   
{
    lnZ = std::log(Zeff);
    //new correction in theta0 formula
    coeffth1 = (1.-8.7780e-2/Zeff)*(0.87+0.03*lnZ);                   
    coeffth2 = (4.0780e-2+1.7315e-4*Zeff)*(0.87+0.03*lnZ);              
    // tail parameters
    G4double lnZ1 = std::log(Zeff+1.);
    coeffc1  = 2.943-0.197*lnZ1;                  
    coeffc2  = 0.0987-0.0143*lnZ1;                              
    // for single scattering
    Z2 = Zeff*Zeff;
    Z23 = std::exp(2.*lnZ/3.);
    scr1     = scr1ini*Z23;
    scr2     = scr2ini*Z2*ChargeSquare;
                                              
    Zold = Zeff;
}

#endif

