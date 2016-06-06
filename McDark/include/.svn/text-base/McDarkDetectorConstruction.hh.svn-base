//-------------------------------------------------------------------
//
// Created by: Timur Sahin (tcsahin@MIT.EDU)
// Date:       July 15, 2007
// Copyright:  MIT 2007
//
//
//
// $Id: McDarkDetectorConstruction.hh,v 1.12 2010/12/29 13:15:06 ddujmic Exp $

#ifndef McDarkDetectorConstruction_H
#define McDarkDetectorConstruction_H 1


#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class McDarkElectricFieldSetup;
class G4VisAttributes;
class McDarkTpcGas;

class McDarkDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    McDarkDetectorConstruction();
    ~McDarkDetectorConstruction();

    G4VPhysicalVolume* Construct();
    
  private:

    // Detector components
    void ConstructVacuumVessel(G4ThreeVector position);
    void ConstructTPC(G4ThreeVector position);
    void ConstructShield(G4ThreeVector position);
    
    // Materials
    void DefineMaterials();
    G4Material *G10,
        *Steel,
        *Epoxy,
        *CH;

    
    McDarkTpcGas *CF4Gas, *XeGas;
    
    G4Material *Cu,
        *Pb,
        *Air,
        *Concrete,
        *Glass,
	*SiliconWafer; // defined by G4

    
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
    G4LogicalVolume* concreteSide_log;
    G4LogicalVolume* concreteFace_log;
    G4LogicalVolume* caseWall_log;
    G4LogicalVolume* caseFloor_log;
    G4LogicalVolume* cf4Gas_log;
    G4LogicalVolume* leadBlock_log;
    G4LogicalVolume* stilt_log;
    G4LogicalVolume* plastic_sheet_log;
    G4LogicalVolume* plastic_mount_log;
    G4LogicalVolume* copper_sheet_log;
    G4LogicalVolume* copper_mount_log;
    G4LogicalVolume* g10mesh_log;
    G4LogicalVolume* wire_container_log_b;
    G4LogicalVolume* wire_container_log_t;
    G4LogicalVolume* large_wire_log;
    G4LogicalVolume* small_wire_log;
    G4LogicalVolume* cover_log;
    G4LogicalVolume* window_log;    

    // Physical volumes
    G4VPhysicalVolume* experimentalHall_phys;
    G4VPhysicalVolume* expZone_phys;
    G4VPhysicalVolume* concrete_phys_1;
    G4VPhysicalVolume* concrete_phys_2;
    G4VPhysicalVolume* concrete_phys_3;
    G4VPhysicalVolume* concrete_phys_4;
    G4VPhysicalVolume* caseWall_phys;
    G4VPhysicalVolume* caseFloor_phys;
    G4VPhysicalVolume* cf4Gas_phys;
    G4VPhysicalVolume* plastic_sheet_bottom;
    G4VPhysicalVolume* copper_sheet_bottom;
    G4VPhysicalVolume* plastic_mount_top;
    G4VPhysicalVolume* copper_mount_top;
    G4VPhysicalVolume* g10mesh_phys;
    G4VPhysicalVolume* wire_container_phys_b;
    G4VPhysicalVolume* wire_container_phys_t;
    G4VPhysicalVolume* cover_phys;
    G4VPhysicalVolume* window_phys;




    
    // EM Field Classes
    McDarkElectricFieldSetup* fEmFieldSetup;

    void DefineField();

    double lowerGroundY;
    double anodeX;
    double anodeY;  
    double anodeA;
    double anodeS;
    double upperGroundY;
    double cathodeY;

};

#endif
