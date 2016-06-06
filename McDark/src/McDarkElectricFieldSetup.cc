//
// Created by: T. Sahin (tcsahin@MIT.EDU), D. Dujmic (ddujmic@MIT.EDU)
// Date:       Set 15, 2007
// Copyright:  MIT 2007
//
//
// $Id: McDarkElectricFieldSetup.cc,v 1.5 2008/03/11 20:44:22 ddujmic Exp $
// GEANT4 tag $Name:  $
//
//  
//   User Field class implementation.
//

#include "McDarkElectricFieldSetup.hh"
#include "McDarkElectricFieldMessenger.hh"

#include "G4UniformElectricField.hh"
#include "McDarkElectricFieldMWPC.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

McDarkElectricFieldSetup::McDarkElectricFieldSetup()
  : fChordFinder(0), fStepper(0), fIntgrDriver(0)
{
    fEMfield = new  G4UniformElectricField( G4ThreeVector( 0.0, 0.0, 0.0 ) );
    //fEMfield = new McDarkElectricFieldMWPC;
    //((McDarkElectricFieldMWPC*)fEMfield)->setAnodeVoltage(anodeVoltage);
    //((McDarkElectricFieldMWPC*)fEMfield)->setCathodeVoltage(cathodeVoltage);
    
    fFieldMessenger = new McDarkElectricFieldMessenger(this) ;  
    fEquation = new G4EqMagElectricField(fEMfield); 
    fMinStep     = 1*um ; // minimal step 
    fStepperType = 1;        
    fFieldManager = GetGlobalFieldManager();
    UpdateField();   
}



////////////////////////////////////////////////////////////////////////////////

McDarkElectricFieldSetup::~McDarkElectricFieldSetup()
{
  if(fChordFinder) delete fChordFinder;
  if(fStepper)     delete fStepper;
  if(fEquation)    delete fEquation;   
  if(fEMfield)     delete fEMfield;
}

/////////////////////////////////////////////////////////////////////////////
//
// Register this field to 'global' Field Manager and 
// Create Stepper and Chord Finder with predefined type, minstep (resp.)
//

void McDarkElectricFieldSetup::UpdateField()
{
    
    SetStepper();

    // update field in manager
    fFieldManager->SetDetectorField( fEMfield );

    // set integrator
    fIntgrDriver = new G4MagInt_Driver( fMinStep, 
                                        fStepper, 
                                        fStepper->GetNumberOfVariables() );

    // set chord. finder
    if (fChordFinder) delete fChordFinder;
    fChordFinder = new G4ChordFinder(fIntgrDriver);
    fFieldManager->SetChordFinder( fChordFinder );
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

void McDarkElectricFieldSetup::SetStepper()
{
  G4int nvar = 8;

  if(fStepper) delete fStepper;
  
  switch ( fStepperType ) 
  {
    case 0:  
      fStepper = new G4ExplicitEuler( fEquation, nvar ); 
      G4cout<<"G4ExplicitEuler is called"<<G4endl;     
      break;
    case 1:  
      fStepper = new G4ImplicitEuler( fEquation, nvar );      
      G4cout<<"G4ImplicitEuler is called"<<G4endl;     
      break;
    case 2:  
      fStepper = new G4SimpleRunge( fEquation, nvar );        
      G4cout<<"G4SimpleRunge is called"<<G4endl;     
      break;
    case 3:  
      fStepper = new G4SimpleHeum( fEquation, nvar );         
      G4cout<<"G4SimpleHeum is called"<<G4endl;     
      break;
    case 4:  
      fStepper = new G4ClassicalRK4( fEquation, nvar );       
      G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;     
      break;
    case 5:  
      fStepper = new G4CashKarpRKF45( fEquation, nvar );      
      G4cout<<"G4CashKarpRKF45 is called"<<G4endl;     
      break;
    case 6:  
      fStepper = 0; // new G4RKG3_Stepper( fEquation, nvar );       
      G4cout<<"G4RKG3_Stepper is not currently working for Electric Field"<<G4endl;     
      break;
    case 7:  
      fStepper = 0; // new G4HelixExplicitEuler( fEquation ); 
      G4cout<<"G4HelixExplicitEuler is not valid for Electric Field"<<G4endl;     
      break;
    case 8:  
      fStepper = 0; // new G4HelixImplicitEuler( fEquation ); 
      G4cout<<"G4HelixImplicitEuler is not valid for Electric Field"<<G4endl;     
      break;
    case 9:  
      fStepper = 0; // new G4HelixSimpleRunge( fEquation );   
      G4cout<<"G4HelixSimpleRunge is not valid for Electric Field"<<G4endl;     
      break;
    default: fStepper = 0;
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field to fieldValue 
//

void McDarkElectricFieldSetup::SetFieldValue(G4double fieldValue)
{
    G4ThreeVector fieldVector( 0.0, fieldValue, 0.0 );
    
    SetFieldValue( fieldVector );
}

///////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field to fieldVector
//

void McDarkElectricFieldSetup::SetFieldValue(G4ThreeVector fieldVector)
{
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr= GetGlobalFieldManager();
    
  if(fieldVector != G4ThreeVector(0.,0.,0.)) {

      G4cout << "Setting uniform electric field " << fieldVector << G4endl;
      
    if(fEMfield) delete fEMfield;    
    fEMfield = new  G4UniformElectricField(fieldVector);

    fEquation->SetFieldObj(fEMfield);  // must now point to the new field

    UpdateField();
   
  }
  
  else {
      
    // If the new field's value is Zero, then it is best to
    //  insure that it is not used for propagation.
    if(fEMfield) delete fEMfield;
    fEMfield = 0;
    fEquation->SetFieldObj(fEMfield);   // As a double check ...

    fieldMgr->SetDetectorField(fEMfield);
  }
  
}


////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager*  McDarkElectricFieldSetup::GetGlobalFieldManager() {
    
  return G4TransportationManager::GetTransportationManager()->GetFieldManager();

}


