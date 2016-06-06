//
// $Id: McDarkRunActionMessenger.cc,v 1.4 2009/12/30 22:34:49 ddujmic Exp $
// GEANT4 tag $Name:  $
//
// --------------------------------------------------------------

#include "McDarkRunActionMessenger.hh"

#include <sstream>

#include "McDarkRunAction.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

McDarkRunActionMessenger::McDarkRunActionMessenger(McDarkRunAction* run)
: McDarkRun(run)
{
    PlotEventCmd = new G4UIcmdWithABool("/mcdark/plotEvent",this);
    PlotEventCmd->SetGuidance("Flag for saving event for offline display");
    PlotEventCmd->SetGuidance("Default = false");
    PlotEventCmd->SetParameterName("plotEvent", false);
    PlotEventCmd->SetDefaultValue(false);

    PrintTracksCmd = new G4UIcmdWithAnInteger("/mcdark/printTracks",this);
    PrintTracksCmd->SetGuidance("Flag for printing a table with tracks in event");
    PrintTracksCmd->SetGuidance(" <1           no printing ");
    PrintTracksCmd->SetGuidance(" =1 (Default) print first 10 tracks; skip empty events");
    PrintTracksCmd->SetGuidance(" =2           print all tracks");
    PrintTracksCmd->SetParameterName("printTracks", false);
    PrintTracksCmd->SetDefaultValue(1);

    PrintHitsCmd = new G4UIcmdWithABool("/mcdark/printHits",this);
    PrintHitsCmd->SetGuidance("Flag for printing a table with hits in event");
    PrintHitsCmd->SetGuidance("Default = false");
    PrintHitsCmd->SetParameterName("printHits", false);
    PrintHitsCmd->SetDefaultValue(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

McDarkRunActionMessenger::~McDarkRunActionMessenger()
{
    delete PlotEventCmd;
    delete PrintTracksCmd;
    delete PrintHitsCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void McDarkRunActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 

    if(command == PlotEventCmd) {
        G4int vl;
        const char* t = newValue;
        std::istringstream is(t);
        is >> vl;
        McDarkRun->setPlotEvent(vl!=0);
    }

    if(command == PrintTracksCmd) {
        G4int vl;
        const char* t = newValue;
        std::istringstream is(t);
        is >> vl;
        McDarkRun->setPrintTracks(vl);
    }

    if(command == PrintHitsCmd) {
        G4int vl;
        const char* t = newValue;
        std::istringstream is(t);
        is >> vl;
        McDarkRun->setPrintHits(vl!=0);
    }

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





