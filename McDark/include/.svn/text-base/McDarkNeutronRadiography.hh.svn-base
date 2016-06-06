//-------------------------------------------------------------------
//
// Created by: Denis Dujmic (ddujmic@MIT.EDU)
// Copyright:  MIT 2009
//
//
//
// $Id: McDarkNeutronRadiography.hh,v 1.6 2010/01/06 15:34:34 ddujmic Exp $

#ifndef McDarkNeutronRadiography_H
#define McDarkNeutronRadiography_H 1


#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class McDarkElectricFieldSetup;
class G4VisAttributes;
class McDarkTpcSD;
class McDarkTpcDigitizer;
class McDarkTpcGas;

class McDarkNeutronRadiography : public G4VUserDetectorConstruction {
    
  public:

    McDarkNeutronRadiography();
    ~McDarkNeutronRadiography();

    G4VPhysicalVolume* Construct();    
    
  private:

    // Detector components
    void ConstructTarget();
    void ConstructTpc0(G4String name);
    void ConstructTpc1(G4String name);
    void ConstructExperimentalHall();
    
    // Materials
    void DefineMaterials();
    G4Material *G10,
        *Steel,
        *Ricorad,
        *Epoxy,
        *CH; // defined here

    McDarkTpcGas *CF4Gas, *HeGas;
    
    G4Material *Cu,
        *Pb, *Al,
        *Air,
        *Concrete,
        *Glass; // defined by G4

    
    // Material visibility
    void DefineVisibility();
    G4VisAttributes *concrete_logVisAtt,
        *steel_logVisAtt,
        *gas_logVisAtt,
        *lead_logVisAtt,
        *plastic_logVisAtt,
        *copper_logVisAtt,
        *g10_logVisAtt,
        *top_wires_VisAtt,
        *bottom_wires_VisAtt;   

    
    // Logical volumes
    G4LogicalVolume* experimentalHall_log;
    G4LogicalVolume* expZone_log;
 

    // Physical volumes
    G4VPhysicalVolume* experimentalHall_phys;
    G4VPhysicalVolume* expZone_phys;

};

#endif
