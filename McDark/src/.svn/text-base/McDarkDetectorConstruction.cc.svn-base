//-------------------------------------------------------------------
//
// Created by: Timur Sahin (tcsahin@MIT.EDU)
// Date:       July 15, 2007
// Copyright:  MIT 2007
//
//
//
// $Id: McDarkDetectorConstruction.cc,v 1.22 2011/03/24 18:33:08 jplopez Exp $

#include "McDarkDetectorConstruction.hh"
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

// Geometry Manager
#include "G4GeometryManager.hh"

//Visualization
#include "G4VisAttributes.hh"
#include "G4Colour.hh"


#include "G4SDManager.hh"
#include "G4DigiManager.hh"

#include "MaxCamSRIM.hh"
#include "McDarkTpcDigitizer.hh"
#include "McDarkCamera.hh"
#include "McDarkTpcGas.hh"
#include "McDarkTpcSD.hh"
#include "McDarkTpcDigitizer.hh"
#include "McDarkDigiManager.hh"
#include "McDarkWormDigitizer.hh"

McDarkDetectorConstruction::McDarkDetectorConstruction()
 : experimentalHall_log(0),
   expZone_log(0),
   concreteSide_log(0),
   concreteFace_log(0),
   experimentalHall_phys(0),
   expZone_phys(0),
   concrete_phys_1(0),
   concrete_phys_2(0),
   concrete_phys_3(0),
   concrete_phys_4(0)
{
    DefineVisibility();
    DefineMaterials();
}


McDarkDetectorConstruction::~McDarkDetectorConstruction() {}



void McDarkDetectorConstruction::DefineVisibility() {
    
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


  
void McDarkDetectorConstruction::DefineMaterials() {

    // Necessary G4prims
    G4double density, temperature, pressure;
    G4String name;
    G4int natoms, nelem;
    
    G4NistManager* man = G4NistManager::Instance();
    
    // Accessing built-in elements
    G4Element* H = man->FindOrBuildElement("H");	// Hydrogen
    G4Element* He = man->FindOrBuildElement("He");	// Hydrogen
    G4Element* C = man->FindOrBuildElement("C");	// Carbon
    G4Element* F = man->FindOrBuildElement("F");	// Fluorine
    G4Element* O = man->FindOrBuildElement("O");	// Oxygen
    G4Element* Si = man->FindOrBuildElement("Si");// Silicon
    G4Element* Fe = man->FindOrBuildElement("Fe");
    G4Element* Mn = man->FindOrBuildElement("Mn");
    G4Element* P = man->FindOrBuildElement("P");
    G4Element* S = man->FindOrBuildElement("S");
    G4Element* Cr = man->FindOrBuildElement("Cr");
    G4Element* Ni = man->FindOrBuildElement("Ni");
    G4Element* Ti = man->FindOrBuildElement("Ti");
    G4Element* Xe = man->FindOrBuildElement("Xe");
    
    // Plastic, CH
    CH = new G4Material("Plastic", density=1.04*g/cm3, nelem=2);
    CH->AddElement(C, natoms=1);
    CH->AddElement(H, natoms=1);
    
   // Tetrafluoromethane, CF4
    name="CF4Gas";
    temperature = 293. * kelvin;		// 20 degrees C
    G4double pTorr = 600; // Torr
    pressure = pTorr/760. * atmosphere;
    //molar mass = (0.875 * 4 g/mol) + (.125 * 88 g/mol)
    density = MaxCamSRIM::density(pTorr, 14.5, temperature) * g / cm3;
    CF4Gas = new McDarkTpcGas(name, density, 3, kStateGas, temperature, pressure);
    CF4Gas->AddElement(C, natoms=25);
    CF4Gas->AddElement(F, natoms=100);
    CF4Gas->AddElement(He, natoms=875);
    CF4Gas->setWorkFunction( 34.3*eV );
    CF4Gas->setAvalancheGain( 5e4 );
    CF4Gas->setPhotonsPerElectron( 0.1 );
    G4double  Efield = 125*volt/cm;
    CF4Gas->setTDiffusion( 2*0.05*volt/Efield );
    CF4Gas->setLDiffusion( 2*0.05*volt/Efield ); // same for now, fix later
    CF4Gas->setTSpread( 0.2*mm );
    CF4Gas->setLSpread( 0.2*mm );
    CF4Gas->setDriftV( 1e7*cm/s );

    // Xe
    pressure = 1 * atmosphere;
    density = 5.894e-3 * g / cm3; 
    XeGas = new McDarkTpcGas(name, density, 1, kStateGas, temperature, pressure);
    XeGas->AddElement(Xe, natoms=1);
    XeGas->setWorkFunction( 22*eV );

    
    // G10, SiO2C3H3
    name="G10";
    density = 1.700 * g / cm3;
    nelem = 4;
    G10 = new G4Material(name, density, nelem);
    G10->AddElement(Si, natoms=1);
    G10->AddElement(O, natoms=2);
    G10->AddElement(C, natoms=3);
    G10->AddElement(H, natoms=3);
    
    // Stainless Steel (1.4541, Stahlschlüssel)
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
    
    // Epoxide
    name = "Epoxy";	// using simple epoxide C2H4O
    density = 899 * kg / m3;
    nelem = 3;
    G4Material* Epoxy = new G4Material(name, density, nelem);
    Epoxy->AddElement(C, natoms=2);
    Epoxy->AddElement(H, natoms=4);
    Epoxy->AddElement(O, natoms=1);
    
    
    // Remaining materials built from within Geant4
    // G4Material* Al = man->FindOrBuildMaterial("G4_Al");	// Aluminum
    Cu = man->FindOrBuildMaterial("G4_Cu");		// Copper
    Pb = man->FindOrBuildMaterial("G4_Pb");		// Lead
    Air = man->FindOrBuildMaterial("G4_AIR");		// Air
    Concrete = man->FindOrBuildMaterial("G4_CONCRETE"); // Concrete
    Glass = man->FindOrBuildMaterial("G4_Pyrex_Glass"); // Glass

    SiliconWafer = new G4Material("SiliconWafer", 2.33*g/cm3, 1);
    SiliconWafer->AddElement(Si, 1);

}



  
G4VPhysicalVolume* McDarkDetectorConstruction::Construct() {

  // GEOMETRIES======================================================

  // EXPERIMENT HALL (World Volume) ---------------------------------
  // Geometry: Cube of air (for now)
  // Will contain entire experiment zone, also possibly background
  // sources.

  // Dimensions
  // Note: the dimensions represent HALF the length across that axis.
  G4double expHall_x = 2.0 * m;
  G4double expHall_y = 2.0 * m;
  G4double expHall_z = 2.0 * m;

  // Construction
  G4Box* experimentalHall_sol = new G4Box("expHall_box", expHall_x, expHall_y, expHall_z);
  
  experimentalHall_log = new G4LogicalVolume(experimentalHall_sol,
                                             Air,
                                             "expHall_log",
                                             0,
                                             0,
                                             0);

  experimentalHall_phys = new G4PVPlacement(0,
                                            G4ThreeVector(),
                                            experimentalHall_log,
                                            "expHall",
                                            0,
                                            false,
                                            0);

  // Visual Attributes:
  // The experiment hall will remain invisible for now, unless modified
  // to account for some background source.
  experimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  // EXPERIMENT BOX (inside EXPERIMENT HALL) ------------------------
  // Geometry: Cube of Air (for now)
  // Another cube of air, inside of the experiment hall.
  // Approximate size of primary experiment area, meant for easy
  // movement/replication of entire detector.
  
  // Dimensions
  // Note: the dimensions represent HALF the length across that axis.
  G4double expZone_side = 1. * m;
  
  // Construction
  G4Box* expZone_sol = new G4Box("expZone_box",
                                 expZone_side,
                                 expZone_side,
                                 expZone_side);
  
  expZone_log = new G4LogicalVolume(expZone_sol,
                                    Air,
                                    "expZone_log");
  
  expZone_phys = new G4PVPlacement(0,
                                   G4ThreeVector(),
                                   "expZone",
                                   expZone_log,
                                   experimentalHall_phys,
                                   false,
                                   0);

  // Visual Attributes:
  // The experiment zone is aphysical and hence invisible.
  expZone_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  

  
  //ConstructVacuumVessel( G4ThreeVector(0,0,0) );


  //ConstructShield( G4ThreeVector(0,0,0) );



  G4double distanceFromBottom=0*cm;
  
  ConstructTPC( G4ThreeVector(0, distanceFromBottom, 0) );



  
  return experimentalHall_phys;
}






void
McDarkDetectorConstruction::ConstructShield(G4ThreeVector pos) {
  // CONCRETE WALL: LONG (inside EXPERIMENT ZONE) -------------------
  // Geometry: Box of Concrete
  // These walsl are in the YZ plane, separated on the X axis.
  
  // Dimensions
  //G4cout << "================= modified side of concrete wall =================" << G4endl;
  G4double concreteSide_x = (9.1 / 2) * cm;
  G4double concreteSide_y = (99.5 / 2) * cm;
  G4double concreteSide_z = (80. / 2) * cm;
  G4double concreteSide_disp = 30.5 * cm;
  
  // Construction
  // NOTE: Using G4PVPlacement as G4Replica requires volume
  // to be the unique daughter of mother daughter, and we
  // need to put other things into 
  G4Box* concreteSide_sol = new G4Box("conceretSide_box",
                                      concreteSide_x,
                                      concreteSide_y,
                                      concreteSide_z);
  
  concreteSide_log = new G4LogicalVolume(concreteSide_sol,
                                         Concrete,
                                         "concreteSide_log");

  // +x
  concrete_phys_1 = new G4PVPlacement(0,
                                      G4ThreeVector(concreteSide_disp, 0, 0),
                                      "concreteSide",
                                      concreteSide_log,
                                      expZone_phys,
                                      false,
                                      0);
  // -x
  concrete_phys_2 = new G4PVPlacement(0,
                                      G4ThreeVector((-1 * concreteSide_disp), 0, 0),
                                      "concreteSide",
                                      concreteSide_log,
                                      expZone_phys,
                                      false,
                                      1);

  // Visual Attributes
  concreteSide_log->SetVisAttributes(concrete_logVisAtt);

  // CONCRETE WALL: SHORT (inside EXPERIMENT ZONE) ------------------
  // Geometry: Box of Concrete
  // These walsl are in the XY plane, separated on the Z axis.
  
  // Dimensions
  G4double concreteFace_x = (52. / 2) * cm;
  G4double concreteFace_y = (99.5 / 2) * cm;
  G4double concreteFace_z = (9.1 / 2) * cm;
  G4double concreteFace_disp = 35.5 * cm;


  // Construction
  G4Box* concreteFace_sol = new G4Box("concreteSide_face",
                                      concreteFace_x,
                                      concreteFace_y,
                                      concreteFace_z);
  concreteFace_log = new G4LogicalVolume(concreteFace_sol,
                                         Concrete,
                                         "concreteFace_log");

  // +z
  concrete_phys_3 = new G4PVPlacement(0,
                                      G4ThreeVector(0, 0, concreteFace_disp),
                                      "concreteFace",
                                      concreteFace_log,
                                      expZone_phys,
                                      false,
                                      0);
  // -z
  concrete_phys_4 = new G4PVPlacement(0,
                                      G4ThreeVector(0, 0, (-1 * concreteFace_disp)),
                                      "concreteFace",
                                      concreteFace_log,
                                      expZone_phys,
                                      false,
                                      1);

  // Visual Attributes
  concreteFace_log->SetVisAttributes(concrete_logVisAtt);
}





void
McDarkDetectorConstruction::ConstructVacuumVessel(G4ThreeVector pos) {


    // STEEL CASE WALL (inside EXPERIMENT ZONE) -----------------------
    // Geometry: Tube of Aluminum
    // Casing used to house the rest of the interior of the detector.

    // Generate rotation matrix 90 degrees about X-axis
    // and its corresponding inverse matrix.
    G4RotationMatrix* rm = new G4RotationMatrix; 
    rm->rotateX(-90 * deg); 
    G4RotationMatrix* antirm = new G4RotationMatrix;
    antirm->rotateX(90 * deg); 

    // Dimensions
    G4double caseWall_rout = 29.6 * cm;
    G4double caseWall_rin = 30.5 * cm;
    G4double caseWall_h = (76.2 / 2) * cm;
    G4double caseWall_startAng = 0 * deg;
    G4double caseWall_endAng = 360 * deg;
    G4double caseWall_disp = 9 * cm; // Displacement from center of
				      // experiment zone

    // Construction
    G4Tubs* caseWall_sol = new G4Tubs("caseWall", 
                                      caseWall_rin, 
                                      caseWall_rout,
                                      caseWall_h, 
                                      caseWall_startAng, 
                                      caseWall_endAng);

    caseWall_log = new G4LogicalVolume(caseWall_sol, 
                                       Steel, 
                                       "caseWall_log");

    G4ThreeVector caseWall_pos = G4ThreeVector(0, caseWall_disp, 0) + pos;
    caseWall_phys = new G4PVPlacement(rm, 
                                      caseWall_pos, 
                                      "casewall", 
                                      caseWall_log,
                                      expZone_phys,
                                      false,
                                      0);
    

    // Visual Attributes
    caseWall_log->SetVisAttributes(steel_logVisAtt);



    
    // STEEL CASE FLOOR (inside CASE WALL) ----------------------------
    // Geometry: Tube cylinder at base of case wall above.
 
    // Dimensions
    G4double caseFloor_rout = 29.6 * cm;
    G4double caseFloor_rin = 0 * cm;
    G4double caseFloor_h = 0.5 * cm;
    G4double caseFloor_startAng = 0 * deg;
    G4double caseFloor_endAng = 360 * deg;
    G4double caseFloor_disp = (-8 * cm) + caseFloor_h; 

    // Construction
    G4Tubs* caseFloor_sol = new G4Tubs("caseFloor", 
                                       caseFloor_rin,
                                       caseFloor_rout,
                                       caseFloor_h,
                                       caseFloor_startAng, 
                                       caseFloor_endAng);
  
    caseFloor_log = new G4LogicalVolume(caseFloor_sol, 
                                        Steel, 
                                        "casefloor_log");

    caseFloor_phys = new G4PVPlacement(rm, 	// No rotation with respect to wall geometry
                                       G4ThreeVector(0, caseFloor_disp, 0), 
                                       "casefloor", 
                                       caseFloor_log,
                                       expZone_phys, 
                                       false, 
                                       0);

    // Visual Attributes
    caseFloor_log->SetVisAttributes(steel_logVisAtt);


    // CF4 gas
    G4Tubs*  vesselGas_geo = new G4Tubs("vesselGas_geo", 0, caseWall_rin, caseWall_h-caseFloor_h, 0*deg, 360*deg);
    G4LogicalVolume *vesselGas_log = new G4LogicalVolume( vesselGas_geo, CF4Gas, "vesselGas_log");
    new G4PVPlacement(antirm, caseWall_pos, "vesselGas", vesselGas_log, expZone_phys, false, 0);




    
    // STEEL CYLINDER COVER (inside EXPERIMENT ZONE) ------------------
    // Geometry: Boolean combination of three cylinder tubes
    
    // Dimensions
    //   Large disk part
    G4double cover_disk_rin = (15.4 / 2) * cm;
    G4double cover_disk_rout = 23 * cm;
    G4double cover_disk_h = (1.0 / 2) * cm; // need half height
    //   Top extension
    G4double cover_ext_rin = cover_disk_rin;
    G4double cover_ext_rout = cover_ext_rin + (1 * cm);
    G4double cover_ext_h = (2.6 / 2) * cm;
    //   Protruding part
    G4double cover_prod_rin = 7.2 * cm;
    G4double cover_prod_rout = 9 * cm;
    G4double cover_prod_h = (2.9 / 2) * cm;
    //   General
    G4double cover_startAng = 0 * deg;
    G4double cover_endAng = 360 * deg;
    G4double cover_disp = caseWall_disp + caseWall_h + cover_disk_h;
    
    
    // Construction
    G4Tubs* cover_disk_sol = new G4Tubs("cover_disk_sol",
                                        cover_disk_rin,
                                        cover_disk_rout,
                                        cover_disk_h,
                                        cover_startAng,
                                        cover_endAng);
    
    G4Tubs* cover_ext_sol = new G4Tubs("cover_ext_sol",
                                       cover_ext_rin,
                                       cover_ext_rout,
                                       cover_ext_h,
                                       cover_startAng,
                                       cover_endAng);
    
    G4Tubs* cover_prod_sol = new G4Tubs("cover_prod_sol",
                                        cover_prod_rin,
                                        cover_prod_rout,
                                        cover_prod_h,
                                        cover_startAng,
                                        cover_endAng);
    
    G4ThreeVector zTrans(0, 0, (cover_disk_h + cover_ext_h));
    G4RotationMatrix* rot = new G4RotationMatrix;
    
    G4UnionSolid* basic_cover_sol = new G4UnionSolid("basic_cover_sol", 
                                                     cover_disk_sol,
                                                     cover_ext_sol,
                                                     rot,
                                                     zTrans);
    
    G4ThreeVector zTrans2(0, 0, (cover_disk_h + cover_ext_h + cover_prod_h));
    
    G4UnionSolid* cover_solid = new G4UnionSolid("cover_solid",
                                                 basic_cover_sol,
                                                 cover_prod_sol,
                                                 rot,
                                                 zTrans2);
    
    cover_log = new G4LogicalVolume(cover_solid, Steel, "case_cover");

    cover_phys = new G4PVPlacement(antirm, 
                                   G4ThreeVector(0, cover_disp, 0),
                                   "cover",
                                   cover_log,
                                   expZone_phys,
                                   false,
                                   0);    // top cover

    
    // Visual Attributes:
    cover_log->SetVisAttributes(steel_logVisAtt);
    
    // WINDOW (inside COVER) ------------------------------------------
    // Geometry: Cylinder tube of glass
    
    G4double win_rin = 0 * cm;
    G4double win_rout = cover_prod_rin;
    G4double win_h = (0.5 / 2) * cm;
    G4double win_startAng = 0 * deg;
    G4double win_endAng = 360 * deg;
    
    G4Tubs* window_sol = new G4Tubs("window_sol", 
                                    win_rin, 
                                    win_rout, 
                                    win_h, 
                                    win_startAng, 
                                    win_endAng);
    
    window_log = new G4LogicalVolume(window_sol,
                                     Glass,
                                     "window");

    window_phys = new G4PVPlacement(rm,
                                    G4ThreeVector(0, 
                                                  cover_disp
                                                  + cover_ext_h 
                                                  + (2 * cover_prod_h)
                                                  + win_h,
                                                  0
                                                  ),
                                    "window",
                                    window_log,
                                    expZone_phys,
                                    false,
                                    0);
}




void
McDarkDetectorConstruction::ConstructTPC(G4ThreeVector position) {

    
    G4String dname="TPC";

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
    G4RotationMatrix* antirm = new G4RotationMatrix;
    antirm->rotateX( 90 * deg);
    
    
    // TPC: gas + field rings
    G4Tubs*  tpc_geo = new G4Tubs(dname+"_geo", 0, 30.5*cm, (63.5 /2)*cm, 0*deg, 360*deg);//inner active volume
    G4LogicalVolume *tpc_log = new G4LogicalVolume( tpc_geo, CF4Gas, dname+"_log");
    G4ThreeVector tpc_pos(position);
    G4PVPlacement *tpc = new G4PVPlacement(antirm, tpc_pos, dname, tpc_log, expZone_phys, false, 0);

    G4Tubs* ring_sol = new G4Tubs(dname+"_ring_sol", 21.08*cm, 22.61*cm, (3.175/2)*mm, 0*deg, 360*deg); //ring inner and outer radius, and height   
    G4LogicalVolume* ring_log = new G4LogicalVolume(ring_sol, Cu, dname+"_ring_log");
    ring_log->SetVisAttributes(copper_logVisAtt);
    
    //number of rings (cathode is considered a ring; anode is not)
    for (G4int i=-27; i<=27; i++) {
        G4ThreeVector pos = tpc_pos + i*G4ThreeVector( 0, (0.71+0.3175)*cm, 0);//ring spacing
        new G4PVPlacement(antirm, pos, dname+"_ring", ring_log, expZone_phys, false, 0);
    }
    G4cout << "Created TPC detector " << tpc->GetName() << G4endl;
 

    // TPC: gaseous sensitive volume
    McDarkTpcSD* tpc_SD = new McDarkTpcSD(dname+"_SD");
    SDman->AddNewDetector( tpc_SD );
    tpc_log->SetSensitiveDetector( tpc_SD );


    // TPC: digitization module
    McDarkTpcDigitizer *tpc_DM = new McDarkTpcDigitizer( dname+"_DM", (const G4PVPlacement*) tpc );
    McDarkDigiManager::addNewModule(tpc_DM);


    // TPC: camera 0 for imaging
    G4Box*  tpc_cam0_geo = new G4Box("tpc_cam0_geo", 1024*12*um, 1024*12*um, 3*um);
    G4LogicalVolume *tpc_cam0_log = new G4LogicalVolume( tpc_cam0_geo, SiliconWafer, dname+"_cam0_log");
    G4ThreeVector tpc_cam0_loc = tpc_pos + (72.0/2)*cm*G4ThreeVector(0,1,0);  // chip distance from amplification plane (focal length)
    G4PVPlacement* tpc_cam0_phys = new G4PVPlacement(antirm, tpc_cam0_loc, dname+"_cam0", tpc_cam0_log, expZone_phys, false, 0);
    McDarkCamera *tpc_cam0=new McDarkCamera(dname+"_cam0", "KAF1001", 500.0, 7.0, 4, 4);
    tpc_cam0->setViewfield(G4ThreeVector(-14.65*cm, -14.65*cm, 0),
                           G4ThreeVector(29.3*cm, 0, 0),
                           G4ThreeVector(0, 29.3*cm, 0),
			   1e-4);
    tpc_cam0->setActiveVolume(tpc_cam0_phys);
    tpc_DM->addCamera( tpc_cam0 );


    // Worms: camera 0 as sensitive volume
    McDarkTpcSD* worm0_SD = new McDarkTpcSD(dname+"_worm0_SD");
    SDman->AddNewDetector( worm0_SD );
    tpc_cam0_log->SetSensitiveDetector( worm0_SD );
    

    // Worms: camera 0 worm digitizer
    McDarkWormDigitizer *worm0_DM = new McDarkWormDigitizer( dname+"_worm0_DM", tpc_cam0 );
    McDarkDigiManager::addNewModule(worm0_DM);

    


    

}



