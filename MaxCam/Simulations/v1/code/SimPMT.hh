/** to do:
    - add gain curve -- gain vs. voltage 
    (exponential by default and user can set the normalization at some voltage)
    - 
 */

#ifndef __SIMPMT__
#define __SIMPMT__

// ROOT Includes
#include "TObject.h"
#include "TObjArray.h"
#include "TString.h"
//#include "TMath.h"
//#include "TRandom3.h"
//#include "TH1.h"

// Viper Includes
#include "SimScope.hh"

// MaxCam Includes
#include "../../../DmtpcPMTQE.hh"
#include "../../../ScopeWaveformData.hh"

class SimPMT : public TObject {

public: 

  SimPMT();
  SimPMT(double x, double y, double diam, double gain, int pmtnum);
  //SimPMT(const SimPMT &other);
  virtual ~SimPMT();

  void setSerialNumber(TString ser){fSerialNumber = ser;}
  void setPMTNumber(Int_t num){fPMTNumber = num;}
  void setDiameter(Double_t diam){fDiameter=diam;}
  void setDistance(Double_t dist) {fDistance=dist;}
  void setPosition(Double_t x, Double_t y);
  void setGain(Double_t gain){fGain = gain;}

  TString  getSerialNumber(){return fSerialNumber;}
  Int_t    getPMTNumber(){return fPMTNumber;}
  Double_t getPositionX(){return fX;}
  Double_t getPositionY(){return fY;}
  Double_t getDiameter() {return fDiameter;}
  Double_t getDistance() {return fDistance;}
  Double_t getGain(){return fGain;}

  DmtpcPMTQE* getPMTQE() {return fPMTQE;}
  void  setPMTQE(DmtpcPMTQE* qe) { fPMTQE = qe; }
  void  setQE(TString filename="", Double_t val=-1.0, TRandom3 *rnd=0);

  TObjArray* wfs() {return fWfs;}
  ScopeWaveformData* wf(Int_t itrig);
  void resetWfs();
  void setWfs() { setWfs(fScope, fScopeChan);}
  void setWfs(SimScope* scope, Int_t chan);

  void assignScope(SimScope* scope, Int_t chan) {fScope=scope; fScopeChan=chan;}
  void      setScope(SimScope* scope) {fScope=scope;}
  SimScope* getScope() {return fScope;}
  void      setScopeChannel(Int_t ch) {fScopeChan=ch;}
  Int_t     getScopeChannel() {return fScopeChan;}

private:

  TString     fSerialNumber; //PMT Serial Number
  Int_t       fPMTNumber;    //PMT Number
  SimScope*   fScope;        //Pointer to scope that this PMT plugs into
  Int_t       fScopeChan;    //Scope channel number that this PMT plugs into
  Double_t    fX;            //[mm] x-Pos of sub-PMT pt of PMT ctr on anode
  Double_t    fY;            //[mm] y-Pos of sub-PMT pt of PMT ctr on anode
  Double_t    fDiameter;     //[mm] Diameter of PMT in mm
  Double_t    fDistance;     //[mm] Perp. dist. from PMT face to sub-PMT pt on the anode. Assumes that PMT face is parallel to anode
  Double_t    fGain;         //Total gain in e- per photon
  DmtpcPMTQE* fPMTQE;        //PMT Quantum Efficiency as a function of wavelength
  Double_t    fWlMin;        //Minimum wavelength for PMT characterization (e.g. QE) [nm]
  Double_t    fWlMax;        //Maximum wavelength for PMT characterization [nm]
  Double_t    fWlStep;       //Resolution in wavelength space for PMT characterization [nm]
  TObjArray*  fWfs;         //Array of waveforms in an event, one per trigger

  ClassDef(SimPMT,1)
};

#endif
