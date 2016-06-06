//
// $Id: McDarkPrimaryGeneratorActionMessenger.cc,v 1.9 2011/01/03 14:14:24 ddujmic Exp $
// GEANT4 tag $Name:  $
//
// --------------------------------------------------------------

#include "McDarkPrimaryGeneratorActionMessenger.hh"
#include "McDarkPrimaryGeneratorAction.hh"

#include <sstream>
#include "assert.h"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "McDarkSpergelDistribution.hh"

McDarkPrimaryGeneratorActionMessenger::McDarkPrimaryGeneratorActionMessenger(McDarkPrimaryGeneratorAction* pga)
: McDarkPrimaryGenerator(pga) {
    
    generatorDirectory = new G4UIdirectory("/generator/");
    generatorDirectory->SetGuidance("Generator control commands");

    spergelDirectory = new G4UIdirectory("/generator/spergel/");
    spergelDirectory->SetGuidance("Commands for Spergel distribution");

    isotropicDirectory = new G4UIdirectory("/generator/isotropic/");
    isotropicDirectory->SetGuidance("Commands for isotropic distribution");
    
    wimpMassCmd = new G4UIcmdWithADouble("/generator/spergel/wimpMass",this);
    wimpMassCmd->SetGuidance("WIMP mass in GeV");
    wimpMassCmd->SetParameterName("wimpMassValue", false);
    wimpMassCmd->SetDefaultValue(100);
    McDarkPrimaryGenerator->getSpergelFunction()->setDarkMatterMass( 100 );

    wimpDensityCmd = new G4UIcmdWithADouble("/generator/spergel/wimpDensity",this);
    wimpDensityCmd->SetGuidance("WIMP density in GeV/c2/cm3");
    wimpDensityCmd->SetParameterName("wimpDensityValue", false);
    wimpDensityCmd->SetDefaultValue(0.3);
    McDarkPrimaryGenerator->getSpergelFunction()->setDarkMatterDensity( 0.3 );
    
    mwVelocityCmd = new G4UIcmdWithADouble("/generator/spergel/mwVelocity",this);
    mwVelocityCmd->SetGuidance("Milky Way veolocity in m/s");
    mwVelocityCmd->SetParameterName("mwVelocityValue", false);
    mwVelocityCmd->SetDefaultValue(230e3);
    McDarkPrimaryGenerator->getSpergelFunction()->setDarkMatterAvgVelocity( 230e3 );

    earthVelocityCmd = new G4UIcmdWithADouble("/generator/spergel/earthVelocity",this);
    earthVelocityCmd->SetGuidance("Earth veolocity around Sun in m/s");
    earthVelocityCmd->SetParameterName("earthVelocityValue", false);
    earthVelocityCmd->SetDefaultValue(244e3);
    McDarkPrimaryGenerator->getSpergelFunction()->setEarthVelocityAroundSun( 244e3 );
    
    earthDVelocityCmd = new G4UIcmdWithADouble("/generator/spergel/earthDVelocity",this);
    earthDVelocityCmd->SetGuidance("Change in Earth veolocity around Sun in m/s");
    earthDVelocityCmd->SetParameterName("earthDVelocityValue", false);
    earthDVelocityCmd->SetDefaultValue(15e3);
    McDarkPrimaryGenerator->getSpergelFunction()->setEarthVelocityVariation( 15e3 );

    yearFractionCmd = new G4UIcmdWithADouble("/generator/spergel/yearFraction",this);
    yearFractionCmd->SetGuidance("Fraction of the year");
    yearFractionCmd->SetParameterName("yearFractionValue", false);
    yearFractionCmd->SetDefaultValue(15e3);
    McDarkPrimaryGenerator->getSpergelFunction()->setEarthVelocityVariation( 15e3 );

    distributionCmd = new G4UIcmdWithAString("/generator/distribution",this);
    distributionCmd->SetGuidance("Energy, angle distribution of recoil nuclei.");
    distributionCmd->SetParameterName("distribution",true);
    distributionCmd->SetDefaultValue("gun");

    minEnergyCmd = new G4UIcmdWithADoubleAndUnit("/generator/minEnergy", this);
    minEnergyCmd->SetGuidance("Minimum energy of flat distribution.");
    //minEnergyCmd->SetParameterName("minEnergyVal");
    minEnergyCmd->SetDefaultValue(0);

    maxEnergyCmd = new G4UIcmdWithADoubleAndUnit("/generator/maxEnergy", this);
    maxEnergyCmd->SetGuidance("Maximum energy of flat distribution.");
    //maxEnergyCmd->SetParameterName("maxEnergyVal");
    maxEnergyCmd->SetDefaultValue(250 * keV);

    minXposCmd = new G4UIcmdWithADoubleAndUnit("/generator/minX", this);
    minXposCmd->SetGuidance("Minimum x position for spatially flat distribution");
    minXposCmd->SetDefaultValue(0 * cm);

    maxXposCmd = new G4UIcmdWithADoubleAndUnit("/generator/maxX",this);
    maxXposCmd->SetGuidance("Maximum x position for spatially flat distribution");
    maxXposCmd->SetDefaultValue(0 * cm);

    minYposCmd = new G4UIcmdWithADoubleAndUnit("/generator/minY", this);
    minYposCmd->SetGuidance("Minimum y position for spatially flat distribution");
    minYposCmd->SetDefaultValue(0 * cm);

    maxYposCmd = new G4UIcmdWithADoubleAndUnit("/generator/maxY", this);
    maxYposCmd->SetGuidance("Maximum y position for spatially flat distribution");
    maxXposCmd->SetDefaultValue(0 * cm);

    minZposCmd = new G4UIcmdWithADoubleAndUnit("/generator/minZ", this);
    minZposCmd->SetGuidance("Minimum z position for spatially flat distribution");
    minZposCmd->SetDefaultValue(0 * cm);

    maxZposCmd = new G4UIcmdWithADoubleAndUnit("/generator/maxZ", this);
    maxZposCmd->SetGuidance("Maximum z position for spatially flat distribution");
    maxZposCmd->SetDefaultValue(0 * cm);
}

McDarkPrimaryGeneratorActionMessenger::~McDarkPrimaryGeneratorActionMessenger() {
    delete wimpMassCmd;
    delete wimpDensityCmd;
    delete mwVelocityCmd;
    delete earthVelocityCmd;
    delete earthDVelocityCmd;
    delete yearFractionCmd;
    delete distributionCmd;

    delete minEnergyCmd;
    delete maxEnergyCmd;
    delete minXposCmd;
    delete maxXposCmd;
    delete minYposCmd;
    delete maxYposCmd;
    delete minZposCmd;
    delete maxZposCmd;
}


void McDarkPrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue) {
    
    if(command == wimpMassCmd) {
        McDarkPrimaryGenerator->getSpergelFunction()->setDarkMatterMass( StoD(newValue) );
    }
    else if(command == wimpDensityCmd) {
        McDarkPrimaryGenerator->getSpergelFunction()->setDarkMatterDensity( StoD(newValue) );
    }
    else if(command == mwVelocityCmd) {
        McDarkPrimaryGenerator->getSpergelFunction()->setDarkMatterAvgVelocity( StoD(newValue) );
    }
    else if(command == earthVelocityCmd) {
        McDarkPrimaryGenerator->getSpergelFunction()->setEarthVelocityAroundSun( StoD(newValue) );
    }
    else if(command == earthDVelocityCmd) {
        McDarkPrimaryGenerator->getSpergelFunction()->setEarthVelocityVariation( StoD(newValue) );
    }
    else if(command == yearFractionCmd) {
        McDarkPrimaryGenerator->getSpergelFunction()->seTimeOfYear( StoD(newValue) );
    }
    else if(command == distributionCmd) {
        newValue.toLower();
        if (newValue=="gun") McDarkPrimaryGenerator->distributionType=McDarkPrimaryGeneratorAction::Gun;
        else if (newValue=="spergel") McDarkPrimaryGenerator->distributionType=McDarkPrimaryGeneratorAction::Spergel;
	else if (newValue=="isotropic") McDarkPrimaryGenerator->distributionType=McDarkPrimaryGeneratorAction::Isotropic;
	else if (newValue=="cf252") McDarkPrimaryGenerator->distributionType=McDarkPrimaryGeneratorAction::Cf252;
	else if (newValue=="dt") McDarkPrimaryGenerator->distributionType=McDarkPrimaryGeneratorAction::DT;
	else if (newValue=="dd") McDarkPrimaryGenerator->distributionType=McDarkPrimaryGeneratorAction::DD;
	else if (newValue=="co57") McDarkPrimaryGenerator->distributionType=McDarkPrimaryGeneratorAction::Co57;
	else assert(0);
    }
    else if(command == minEnergyCmd) {
	McDarkPrimaryGenerator->SetMinEnergy(maxEnergyCmd->GetNewDoubleValue(newValue));
    }
    else if(command == maxEnergyCmd) {
	McDarkPrimaryGenerator->SetMaxEnergy(maxEnergyCmd->GetNewDoubleValue(newValue));
    }
    else if(command == minXposCmd) {
	McDarkPrimaryGenerator->SetMinXPos(minXposCmd->GetNewDoubleValue(newValue));
    }
    else if(command == maxXposCmd) {
	McDarkPrimaryGenerator->SetMaxXPos(maxXposCmd->GetNewDoubleValue(newValue));
    }
    else if(command == minYposCmd) {
	McDarkPrimaryGenerator->SetMinYPos(minYposCmd->GetNewDoubleValue(newValue));
    }
    else if(command == maxYposCmd) {
	McDarkPrimaryGenerator->SetMaxYPos(maxYposCmd->GetNewDoubleValue(newValue));
    }
    else if(command == minZposCmd) {
	McDarkPrimaryGenerator->SetMinZPos(minZposCmd->GetNewDoubleValue(newValue));
    }
    else if(command == maxZposCmd) {
	McDarkPrimaryGenerator->SetMaxZPos(maxZposCmd->GetNewDoubleValue(newValue));
    }

    
}








