//
// Created by: T. Sahin (tcsahin@MIT.EDU), D. Dujmic (ddujmic@MIT.EDU)
// Date:       Set 15, 2007
// Copyright:  MIT 2007
//
//
// $Id: McDarkElectricFieldSetup.hh,v 1.3 2008/03/11 20:44:22 ddujmic Exp $
// GEANT4 tag $Name:  $
//
//    A class for control of the Electric Field of the detector.
//
//


#ifndef McDarkElectricFieldSetup_H
#define McDarkElectricFieldSetup_H

#include "McDarkElectricFieldMWPC.hh"
#include "G4ThreeVector.hh"

class G4FieldManager;
class G4ChordFinder;
class G4EquationOfMotion;
class G4Mag_EqRhs;
class G4EqMagElectricField;
class G4MagIntegratorStepper;
class G4MagInt_Driver; 
class McDarkElectricFieldMessenger;
class G4ElectricField;

class McDarkElectricFieldSetup {
    
public:

  McDarkElectricFieldSetup() ;

 ~McDarkElectricFieldSetup() ;  
      
  void SetStepperType( G4int i) { fStepperType = i ; }

  void SetStepper();

  void SetMinStep(G4double s) { fMinStep = s ; }
  void SetAnodeVoltage(G4double s) { anodeVoltage = s ; }
  void SetCathodeVoltage(G4double s) { cathodeVoltage = s ; }

  void UpdateField();

    McDarkElectricFieldMWPC* getFieldMWPC() { return (McDarkElectricFieldMWPC*)fEMfield; }
    G4ElectricField* getField() { return fEMfield; }
    
  void SetFieldValue(G4ThreeVector fieldVector);
  void SetFieldValue(G4double      fieldValue);
  G4ThreeVector GetConstantFieldValue();

protected:

      // Find the global Field Manager
  G4FieldManager*         GetGlobalFieldManager() ;

private:
    
    G4FieldManager*         fFieldManager ;

    G4ChordFinder*          fChordFinder ;

    G4EqMagElectricField*   fEquation ;

    //McDarkElectricFieldMWPC*        fEMfield;
    G4ElectricField*        fEMfield;
 
    G4ThreeVector           fElFieldValue ; 

    G4MagIntegratorStepper* fStepper ;
    G4MagInt_Driver*        fIntgrDriver;

    G4int                   fStepperType ;

    G4double                fMinStep ;

    G4double  anodeVoltage;
    
    G4double  cathodeVoltage;
    
    McDarkElectricFieldMessenger*      fFieldMessenger;

};

#endif
