//
//
// $Id: McDarkRunAction.hh,v 1.13 2009/12/30 22:34:49 ddujmic Exp $
// GEANT4 tag $Name:  $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef McDarkRunAction_h
#define McDarkRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TCollection.h"

#include "McDarkCameraManager.hh"

class McDarkRunActionMessenger;
class DmtpcDataset;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

class McDarkRunAction : public G4UserRunAction {
    
  public:
    McDarkRunAction();
   ~McDarkRunAction();

  public:

    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

 
    McDarkCameraManager* getCamMan() { return camMan; }


    DmtpcDataset *data() { return _data; }

    
    void setPrintTracks(G4int p) { _printTracks=p; }
    G4int getPrintTracks() { return _printTracks; }
    void setPrintHits(G4bool p) { _printHits=p; }
    G4bool getPrintHits() { return _printHits; }
    void setPlotEvent(G4bool p) { _plotEvent=p; }
    G4bool getPlotEvent() { return _plotEvent; }
    
    
private:
    
    McDarkRunActionMessenger* runMessenger;
    McDarkCameraManager* camMan;

    DmtpcDataset *_data;


    // display and saving
    G4int _printTracks;
    G4bool _printHits;
    G4bool _plotEvent;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif





