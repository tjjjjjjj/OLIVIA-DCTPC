//
// Created by: T. Sahin (tcsahin@MIT.EDU), D. Dujmic (ddujmic@MIT.EDU)
// Date:       Set 15, 2007
// Copyright:  MIT 2007
//
//
// $Id: McDarkElectricFieldMessenger.cc,v 1.1 2007/10/17 11:51:17 ddujmic Exp $
// GEANT4 tag $Name:  $
//
// 

#include "McDarkElectricFieldMessenger.hh"
#include "McDarkElectricFieldSetup.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//////////////////////////////////////////////////////////////////////////////

McDarkElectricFieldMessenger::McDarkElectricFieldMessenger(McDarkElectricFieldSetup* pEMfield)
  :fEFieldSetup(pEMfield)
{

    
  mcDarkdetDir = new G4UIdirectory("/field/");
  mcDarkdetDir->SetGuidance("McDark field control.");

  stepperCmd = new G4UIcmdWithAnInteger("/field/setStepperType",this);
  stepperCmd->SetGuidance("Select stepper type for electric field");
  stepperCmd->SetParameterName("choice",true);
  stepperCmd->SetDefaultValue(4);
  stepperCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  updateCmd = new G4UIcmdWithoutParameter("/field/update",this);
  updateCmd->SetGuidance("Update field configuration in the detector.");
  updateCmd->SetGuidance("This command MUST be applied before \"beamOn\". ");
  updateCmd->AvailableForStates(G4State_Idle);
      
  elecFieldCmd = new G4UIcmdWithADoubleAndUnit("/field/setUniformEy",this);  
  elecFieldCmd->SetGuidance("Define uniform electric field.");
  elecFieldCmd->SetParameterName("Ey",false,false);
  elecFieldCmd->SetDefaultUnit("volt");
  elecFieldCmd->AvailableForStates(G4State_Idle); 
 
  minStepCmd = new G4UIcmdWithADoubleAndUnit("/field/setMinStep",this);  
  minStepCmd->SetGuidance("Define minimal step");
  minStepCmd->SetGuidance("Magnetic field will be in Z direction.");
  minStepCmd->SetParameterName("minStep",false,false);
  minStepCmd->SetDefaultUnit("mm");
  minStepCmd->AvailableForStates(G4State_Idle);  
       
//   dielectricCmd = new G4UIcmdWithAString("/field/dielectric",this);
//   dielectricCmd->SetGuidance("Select Material of the dielectric.");
//   dielectricCmd->SetParameterName("dielectric",true);
//   dielectricCmd->SetDefaultValue("Xe");
//   dielectricCmd->AvailableForStates(G4State_Idle);

  anodeVoltageCmd = new G4UIcmdWithADoubleAndUnit("/field/anodeVoltage",this);  
  anodeVoltageCmd->SetGuidance("Anode voltage");
  anodeVoltageCmd->SetParameterName("anodeVoltage",false,false);
  anodeVoltageCmd->SetDefaultUnit("volt");
  anodeVoltageCmd->AvailableForStates(G4State_Idle);  

  cathodeVoltageCmd = new G4UIcmdWithADoubleAndUnit("/field/cathodeVoltage",this);  
  cathodeVoltageCmd->SetGuidance("Cathode voltage");
  cathodeVoltageCmd->SetParameterName("cathodeVoltage",false,false);
  cathodeVoltageCmd->SetDefaultUnit("volt");
  cathodeVoltageCmd->AvailableForStates(G4State_Idle);  
}

///////////////////////////////////////////////////////////////////////////////

McDarkElectricFieldMessenger::~McDarkElectricFieldMessenger()
{
  delete mcDarkdetDir;
  delete stepperCmd;
  delete updateCmd;  
  delete elecFieldCmd;
  delete minStepCmd;
  //delete dielectricCmd;
  delete anodeVoltageCmd;
  delete cathodeVoltageCmd;
}

////////////////////////////////////////////////////////////////////////////
//
//

void McDarkElectricFieldMessenger::SetNewValue( G4UIcommand* command, G4String newValue)
{
    
  if( command == stepperCmd ) { 
    fEFieldSetup->SetStepperType(stepperCmd->GetNewIntValue(newValue));
  }
  
  else if( command == updateCmd ) { 
    fEFieldSetup->UpdateField(); 
  }
  
  else if( command == elecFieldCmd ) { 
    fEFieldSetup->SetFieldValue(elecFieldCmd->GetNewDoubleValue(newValue));
  }
  
  else if( command == minStepCmd ) { 
    fEFieldSetup->SetMinStep( minStepCmd->GetNewDoubleValue(newValue));
  }

  else if( command == anodeVoltageCmd ) { 
    fEFieldSetup->SetAnodeVoltage( anodeVoltageCmd->GetNewDoubleValue(newValue));
  }

  else if( command == cathodeVoltageCmd ) { 
    fEFieldSetup->SetCathodeVoltage( cathodeVoltageCmd->GetNewDoubleValue(newValue));
  }
  
}

//
//
/////////////////////////////////////////////////////////////////////////
