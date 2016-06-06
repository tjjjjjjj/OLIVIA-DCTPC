//-------------------------------------------------------------------
//
// Created by: Denis Dujmic (dujmic@MIT.EDU)
// Copyright:  MIT 2009
//
//
//
// $Id: McDarkNeutronRadiography.cc,v 1.9 2010/01/06 15:34:34 ddujmic Exp $

#include "McDarkNeutronRadiography.hh"
#include "globals.hh"

// Materials
#include "G4Material.hh"
#include "G4NistManager.hh"

// Solids
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

// Vector manipulation
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

// Fields
#include "McDarkElectricFieldSetup.hh"
#include "McDarkElectricFieldMWPC.hh"
#include "McDarkTpcGas.hh"

// Geometry Manager
#include "G4GeometryManager.hh"

//Visualization
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// Sensitive volume and hits
#include "McDarkTpcSD.hh"
#include "G4SDManager.hh"
#include "MaxCamSRIM.hh"

// Digitizer 
#include "McDarkTpcDigitizer.hh"
#include "G4DigiManager.hh"
#include "McDarkCamera.hh"




McDarkNeutronRadiography::McDarkNeutronRadiography()
 : experimentalHall_log(0), expZone_log(0),
   experimentalHall_phys(0), expZone_phys(0) {

    DefineVisibility();
    DefineMaterials();
    
}


McDarkNeutronRadiography::~McDarkNeutronRadiography() {}



void McDarkNeutronRadiography::DefineVisibility() {
    
    // VISIBILITY (COLORS) GIVE VALUES IN RGB
    concrete_logVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));		//   grey
    steel_logVisAtt = new G4VisAttributes(G4Colour(.1,.5,.5));	//   black
    gas_logVisAtt = new G4VisAttributes(G4Colour(0.0,0.3,0.7));		//   bluish
    lead_logVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));		//   blue
    plastic_logVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));		//   grey
    plastic_logVisAtt->SetForceSolid(true);
    plastic_logVisAtt->SetForceSolid(true);
    copper_logVisAtt = new G4VisAttributes(G4Colour(1.0, 0.7, 0.0));	//   orange
    copper_logVisAtt->SetForceSolid(true);
    g10_logVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));	//   green
    top_wires_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));	// yellow
    top_wires_VisAtt->SetForceSolid(true);
    bottom_wires_VisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));	// cyan
    bottom_wires_VisAtt->SetForceSolid(true);
}


  
void McDarkNeutronRadiography::DefineMaterials() {

    // Necessary G4prims
    G4double density, temperature, pressure;
    G4String name;
    G4int natoms, nelem;
    
    G4NistManager* man = G4NistManager::Instance();
    
    // Accessing built-in elements
    G4Element* H = man->FindOrBuildElement("H");	// Hydrogen
    G4Element* He = man->FindOrBuildElement("He");	// Helium
    G4Element* B = man->FindOrBuildElement("B");	// Boron
    G4Element* C = man->FindOrBuildElement("C");	// Carbon
    G4Element* F = man->FindOrBuildElement("F");	// Fluorine
    G4Element* O = man->FindOrBuildElement("O");	// Oxygen
    G4Element* Si = man->FindOrBuildElement("Si");      // Silicon
    G4Element* Fe = man->FindOrBuildElement("Fe");
    G4Element* Mn = man->FindOrBuildElement("Mn");
    G4Element* P = man->FindOrBuildElement("P");
    G4Element* S = man->FindOrBuildElement("S");
    G4Element* Cr = man->FindOrBuildElement("Cr");
    G4Element* Ni = man->FindOrBuildElement("Ni");
    G4Element* Ti = man->FindOrBuildElement("Ti");

    
    // Plastic, CH
    CH = new G4Material("Plastic", density=1.04*g/cm3, nelem=2);
    CH->AddElement(C, natoms=1);
    CH->AddElement(H, natoms=1);

    
    // Tetrafluoromethane, CF4
    name="CF4Gas";
    temperature = 293. * kelvin;		// 20 degrees C
    G4double pTorr = 75; // Torr
    pressure = pTorr/760. * atmosphere;
    density = MaxCamSRIM::density(pTorr, 88, temperature) * g / cm3;
    CF4Gas = new McDarkTpcGas(name, density, 2, kStateGas, temperature, pressure);
    CF4Gas->AddElement(C, natoms=1);
    CF4Gas->AddElement(F, natoms=4);
    CF4Gas->setWorkFunction( 34.3*eV );
    CF4Gas->setAvalancheGain( 5e4 );
    CF4Gas->setPhotonsPerElectron( 0.1 );
    G4double  Efield = 280*volt/cm;
    CF4Gas->setTDiffusion( 2*0.05*volt/Efield );
    CF4Gas->setLDiffusion( 2*0.05*volt/Efield ); // same for now, fix later
    CF4Gas->setTSpread( 0.2*mm );
    CF4Gas->setLSpread( 0.2*mm );
    CF4Gas->setDriftV( 1e7*cm/s );

   
    

    // He-4 gas
    name="HeGas";
    pTorr = 700.; // Torr
    pressure = pTorr/760. * atmosphere;
    density = MaxCamSRIM::density(pTorr, 4, 293) * g / cm3;
    HeGas = new McDarkTpcGas(name, density, 1, kStateGas, temperature, pressure);
    HeGas->AddElement(He, natoms=1);
    HeGas->setWorkFunction( 34.3*eV );
    HeGas->setAvalancheGain( 5e4 );
    HeGas->setPhotonsPerElectron( 0.1 );
    HeGas->setTDiffusion( 2*0.05*volt/Efield );
    HeGas->setLDiffusion( 2*0.05*volt/Efield ); // same for now, fix later
    HeGas->setTSpread( 0.2*mm );
    HeGas->setLSpread( 0.2*mm );
    HeGas->setDriftV( 1e7*cm/s );

    
    
    // G10, SiO2C3H3
    name="G10";
    density = 1.700 * g / cm3;
    nelem = 4;
    G10 = new G4Material(name, density, nelem);
    G10->AddElement(Si, natoms=1);
    G10->AddElement(O, natoms=2);
    G10->AddElement(C, natoms=3);
    G10->AddElement(H, natoms=3);

    
    
    // Stainless Steel (1.4541, Stahlschlaessel)
    name="Stainless steel";
    density = 7.9*g/cm3;
    nelem = 9;
    Steel = new G4Material(name, density, nelem);
    Steel->AddElement(Fe, 0.68325);
    Steel->AddElement(C,  0.001);
    Steel->AddElement(Si, 0.01);
    Steel->AddElement(Mn, 0.02);
    Steel->AddElement(P,  0.00045);
    Steel->AddElement(S,  0.0003);
    Steel->AddElement(Cr, 0.18);
    Steel->AddElement(Ni, 0.10);
    Steel->AddElement(Ti, 0.005);


    // Ricorad
    Ricorad = new G4Material("Ricorad", 0.945*g/cm3 , 3);
    Ricorad->AddElement(H, 0.11);
    Ricorad->AddElement(B, 0.02);
    Ricorad->AddElement(C, 0.87);

    
    // Epoxide
    name = "Epoxy";	// using simple epoxide C2H4O
    density = 899 * kg / m3;
    nelem = 3;
    G4Material* Epoxy = new G4Material(name, density, nelem);
    Epoxy->AddElement(C, natoms=2);
    Epoxy->AddElement(H, natoms=4);
    Epoxy->AddElement(O, natoms=1);
    
    
    // Remaining materials built from within Geant4
    Cu = man->FindOrBuildMaterial("G4_Cu");		// Copper
    Al = man->FindOrBuildMaterial("G4_Al");		// Aluminum
    Pb = man->FindOrBuildMaterial("G4_Pb");		// Lead
    Air = man->FindOrBuildMaterial("G4_AIR");		// Air
    Concrete = man->FindOrBuildMaterial("G4_CONCRETE"); // Concrete
    Glass = man->FindOrBuildMaterial("G4_Pyrex_Glass"); // Glass
}



G4VPhysicalVolume*
McDarkNeutronRadiography::Construct() {
  
  ConstructExperimentalHall();
  //ConstructTarget();
  ConstructTpc0("det0");
  ConstructTpc1("det1");

  
  return experimentalHall_phys;
}




void  
McDarkNeutronRadiography::ConstructExperimentalHall() {

  G4double expHall_x = 2.0 * m;
  G4double expHall_y = 2.0 * m;
  G4double expHall_z = 2.0 * m;

  
  // Construct hall
  G4Box* experimentalHall_sol = new G4Box("expHall_box", expHall_x, expHall_y, expHall_z);
  
  experimentalHall_log = new G4LogicalVolume(experimentalHall_sol, Air, "expHall_log");

  experimentalHall_phys = new G4PVPlacement(0,
                                            G4ThreeVector(),
                                            experimentalHall_log,
                                            "expHall",
                                            0,
                                            false,
                                            0);
  experimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());


  
  G4double expZone_side = expHall_x;
  
  // Construct experimental zone
  G4Box* expZone_sol = new G4Box("expZone_box", expZone_side, expZone_side, expZone_side);
  
  expZone_log = new G4LogicalVolume(expZone_sol, Air, "expZone_log");
  
  expZone_phys = new G4PVPlacement(0,
                                   G4ThreeVector(),
                                   "expZone",
                                   expZone_log,
                                   experimentalHall_phys,
                                   false,
                                   0);
  expZone_log->SetVisAttributes(G4VisAttributes::GetInvisible());
}





void
McDarkNeutronRadiography::ConstructTarget() {

    
    G4RotationMatrix* antirm = new G4RotationMatrix;
    antirm->rotateX(0 * deg); 
    
    // A plastic target
    G4Tubs*  target0_geo = new G4Tubs("target0_geo", 0, 5*cm, 5*cm, 0, 360);
    G4LogicalVolume *target0_log = new G4LogicalVolume( target0_geo, CH, "target0_log");
    target0_log->SetVisAttributes( plastic_logVisAtt );
    G4PVPlacement *target0 = new G4PVPlacement(antirm,
                                               G4ThreeVector(0, 0, 0),
                                               "target0",
                                               target0_log,
                                               expZone_phys,
                                               false,
                                               0);
    G4cout<<target0->GetName()<<G4endl;
}




void
McDarkNeutronRadiography::ConstructTpc0(G4String dname) {
    
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4DigiManager * DMman = G4DigiManager::GetDMpointer();
    
    G4RotationMatrix* antirm = new G4RotationMatrix;
    antirm->rotateX(0 * deg);

    
    // CF4 TPC detector
    G4Tubs*  tpc_geo = new G4Tubs(dname+"_geo", 0, 15*cm, 10*cm, 0*deg, 360*deg);
    G4LogicalVolume *tpc_log = new G4LogicalVolume( tpc_geo, CF4Gas, dname+"_log");
    G4ThreeVector tpc_pos(0, 0, 30*cm);
    G4PVPlacement *tpc = new G4PVPlacement(antirm, tpc_pos, dname, tpc_log, expZone_phys, false, 0);

    G4Tubs* ring_sol = new G4Tubs(dname+"_ring_sol", 15*cm, 17*cm, 1*mm, 0*deg, 360*deg);    
    G4LogicalVolume* ring_log = new G4LogicalVolume(ring_sol, Cu, dname+"_ring_log");
    ring_log->SetVisAttributes(copper_logVisAtt);
    for (G4int i=-10; i<11; i++) {
        G4ThreeVector pos = tpc_pos + i*G4ThreeVector( 0, 0, 1*cm);
        new G4PVPlacement(antirm, pos, dname+"_ring", ring_log, expZone_phys, false, 0);
    }
    G4cout << "Created detector " << tpc->GetName() << G4endl;

    
    // create sensitive detector volume
    McDarkTpcSD* tpc_SD = new McDarkTpcSD(dname+"_SD");
    SDman->AddNewDetector( tpc_SD );
    tpc_log->SetSensitiveDetector( tpc_SD );
 

    // create digitization module
    McDarkTpcDigitizer *tpc_DM = new McDarkTpcDigitizer( dname+"_DM", (const G4PVPlacement*) tpc );
    DMman->AddNewModule(tpc_DM);
    McDarkCamera *tpc_cam0=new McDarkCamera(dname+"_cam0", "KAF1001", 0.0, 7.0, 4, 4);
    tpc_cam0->setViewfield(G4ThreeVector(-8*cm, -8*cm, 0),
                            G4ThreeVector(16*cm, 0, 0),
                            G4ThreeVector(0, 16*cm, 0),
                            1e-4);
    tpc_DM->addCamera( tpc_cam0 );
}








void
McDarkNeutronRadiography::ConstructTpc1(G4String dname) {
    // This is a dummy detector for testing
    
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4DigiManager * DMman = G4DigiManager::GetDMpointer();

    G4RotationMatrix* antirm = new G4RotationMatrix;
    antirm->rotateX(0 * deg); 
    
    // He TPC detector
    G4Tubs *tpc_geo = new G4Tubs(dname+"_geo", 0, 40*cm, 10*cm, 0*deg, 360*deg);
    G4LogicalVolume *tpc_log = new G4LogicalVolume( tpc_geo, HeGas, dname+"_log");
    G4ThreeVector tpc_pos(0, 0, 60*cm);
    G4PVPlacement *tpc = new G4PVPlacement(antirm, tpc_pos, dname, tpc_log, expZone_phys, false, 0);
    
    G4Tubs* ring_sol = new G4Tubs(dname+"_ring_sol", 40*cm, 42*cm, 1*mm, 0*deg, 360*deg);    
    G4LogicalVolume* ring_log = new G4LogicalVolume(ring_sol, Cu, dname+"_ring_log");
    ring_log->SetVisAttributes(copper_logVisAtt);
    for (G4int i=-10; i<11; i++) {
        G4ThreeVector pos = tpc_pos + i*G4ThreeVector( 0, 0, 1*cm);
        new G4PVPlacement(antirm, pos, dname+"_ring", ring_log, expZone_phys, false, 0);
    }
    G4cout << "Created detector " << tpc->GetName() << G4endl;

    
    // create sensitive detector volume
    McDarkTpcSD* tpcSD = new McDarkTpcSD(dname+"_SD");
    SDman->AddNewDetector( tpcSD );
    tpc_log->SetSensitiveDetector( tpcSD );
    

    // create digitization modules    
    McDarkTpcDigitizer *tpcDM = new McDarkTpcDigitizer( dname+"_DM", (const G4PVPlacement*) tpc );
    DMman->AddNewModule( tpcDM);
    McDarkCamera *tpc_cam0=new McDarkCamera(dname+"_cam0", "KAF1001", 0.0, 7.0, 4, 4);
    tpc_cam0->setViewfield(G4ThreeVector(-8*cm, -8*cm, 0),
                            G4ThreeVector(16*cm, 0, 0),
                            G4ThreeVector(0, 16*cm, 0),
                            1e-4);
    tpcDM->addCamera( tpc_cam0 );

    

}
