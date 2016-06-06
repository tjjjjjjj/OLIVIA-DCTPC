
//
// $Id: McDarkRunAction.cc,v 1.15 2009/12/30 22:34:49 ddujmic Exp $
// GEANT4 tag $Name:  $
// 

#include "McDarkRunAction.hh"
#include "McDarkRunActionMessenger.hh"

#include "G4Run.hh"
#include "G4DigiManager.hh"
#include "G4SDManager.hh"
#include "TCanvas.h"
#include "TSystem.h"

McDarkRunAction::McDarkRunAction() {
    runMessenger = new McDarkRunActionMessenger(this);    
}


McDarkRunAction::~McDarkRunAction()
{}

#include "DmtpcDataset.hh"

void McDarkRunAction::BeginOfRunAction(const G4Run* aRun) {
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

    _data = new DmtpcDataset;

    
    const char* fileDir = gSystem->Getenv("ROOTOUT");
    if (!fileDir) {
	G4cout << "ROOTOUT file directory not found" << G4endl;
	G4cout << "Setting ROOTOUT directory to ./" << G4endl;
	fileDir = ".";
    }
    TString fullName(fileDir);
    fullName += "/mcdark.root";

    _data->createRootFile( (const char*)fullName.Data(), "recreate" );

    
    G4DigiManager * fDM = G4DigiManager::GetDMpointer();
    G4SDManager* fSD = G4SDManager::GetSDMpointer();
    G4cout << "+----------------------------------------------------+" << G4endl;
    G4cout << "|                                                    |" << G4endl;
    G4cout << "|         R U N  S U M M A R Y                       |" << G4endl;
    G4cout << "|                                                    |" << G4endl;
    G4cout << "|                                                    |" << G4endl;
    G4cout << "|                                                    |" << G4endl;
    G4cout << "| Sensitive detector(s), hit collection(s):          |" << G4endl;
    G4cout << "|                                                    |" << G4endl;
    for (G4int i=0; i<fSD->GetHCtable()->entries(); i++) {
        G4cout << "|        " 
               << fSD->GetHCtable()->GetSDname(i) <<  ",   "
               << fSD->GetHCtable()->GetHCname(i) << G4endl; 
    }
    G4cout << "|                                                    |" << G4endl;
    G4cout << "|                                                    |" << G4endl;
    G4cout << "| Digitization module(s), digi collection(s):        |" << G4endl;
    G4cout << "|                                                    |" << G4endl;
    for (G4int i=0; i<fDM->GetDCtable()->entries(); i++) {
        G4cout << "|        " 
               << fDM->GetDCtable()->GetDMname(i) <<  ",   "
               << fDM->GetDCtable()->GetDCname(i) << G4endl; 
    }
    G4cout << "|                                                    |" << G4endl;
    G4cout << "+----------------------------------------------------+" << G4endl;

    
}


void McDarkRunAction::EndOfRunAction(const G4Run*) {
    G4cout << "McDarkRunAction: save data file " << G4endl;
    data()->file()->ls();
    data()->file()->Write();
    data()->file()->Close();
    G4cout << "Done" << G4endl;
}

