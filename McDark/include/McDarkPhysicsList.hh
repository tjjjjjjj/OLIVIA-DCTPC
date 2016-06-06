//-------------------------------------------------------------------
//
// Created by: T. Sahin (tcsahin@MIT.EDU), D. Dujmic (ddujmic@MIT.EDU)
// Date:       July 15, 2007
// Copyright:  MIT 2007
//
//
//
// $Id: McDarkPhysicsList.hh,v 1.3 2007/08/10 17:03:34 ddujmic Exp $
// GEANT4 tag $Name:  $
//
// 
//
// McDarkPhysicsList
//  Construct/define particles and physics processes
//


#ifndef McDarkPhysicsList_h
#define McDarkPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class McDarkPhysicsList: public G4VUserPhysicsList {
    
public:
    McDarkPhysicsList();
    ~McDarkPhysicsList();
    
protected:
    // Construct particle and physics process
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

    // these methods Construct physics processes and register them
    virtual void ConstructGeneral();
    virtual void ConstructEM();
    virtual void ConstructHad();
    virtual void ConstructOp();

    virtual void AddTransportation();
    
private:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructHadrons();
    void ConstructShortLiveds();

    // verbose
    G4int VerboseLevel;
    G4int OpVerbLevel;
    
};

#endif
