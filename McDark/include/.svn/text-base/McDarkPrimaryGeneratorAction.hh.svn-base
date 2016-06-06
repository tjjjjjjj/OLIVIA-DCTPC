//-------------------------------------------------------------------
//
// Created by: Timur Sahin (tcsahin@MIT.EDU)
// Date:       July 15, 2007
// Copyright:  MIT 2007
//
//
// $Id: McDarkPrimaryGeneratorAction.hh,v 1.14 2010/12/28 14:45:49 tcsahin Exp $
// GEANT4 tag $Name:  $
//

#ifndef McDarkPrimaryGeneratorAction_h
#define McDarkPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "McDarkPrimaryGeneratorActionMessenger.hh"

class G4ParticleGun;
class G4Event;

class McDarkSpergelDistribution;


class McDarkPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    
    McDarkPrimaryGeneratorAction();

    ~McDarkPrimaryGeneratorAction();


  public:

    McDarkSpergelDistribution* getSpergelFunction() { return spergelDist; }
    
    void GeneratePrimaries(G4Event* anEvent);

    enum RequestDistribution {
        Gun,
        Spergel,
	Isotropic,
        Cf252,
        DT,
        DD,
        CosmicMuons,
		Co57
    } distributionType;

    void SetMinEnergy(G4double e) { particle_minEnergy = e; }
    void SetMaxEnergy(G4double e) { particle_maxEnergy = e; }
    void SetMinXPos(G4double x) { particle_minX = x; }
    void SetMaxXPos(G4double x) { particle_maxX = x; }
    void SetMinYPos(G4double y) { particle_minY = y; }
    void SetMaxYPos(G4double y) { particle_maxY = y; }
    void SetMinZPos(G4double z) { particle_minZ = z; }
    void SetMaxZPos(G4double z) { particle_maxZ = z; }
    
    void print();
    
  private:
    G4ParticleGun* particleGun;
    
    McDarkSpergelDistribution* spergelDist;

    G4String GetName() { return "McDarkPrimaryGeneratorAction"; }
    
    McDarkPrimaryGeneratorActionMessenger* theMessenger;


    // For isotropic distribution:
    G4double particle_minEnergy;
    G4double particle_maxEnergy;
    G4double particle_minX;
    G4double particle_maxX;
    G4double particle_minY;
    G4double particle_maxY;
    G4double particle_minZ;
    G4double particle_maxZ;


};

#endif


