#ifndef __SCOPEDATAINFO__
#define __SCOPEDATAINFO__

#include "TObject.h"
#include "TObjArray.h"
#include "TString.h"

class ScopeDataInfo : public TObject {

public: 

  ScopeDataInfo();
  ScopeDataInfo(Int_t num, TString name, TString cd);
  ScopeDataInfo(const ScopeDataInfo &other);
  virtual ~ScopeDataInfo();

  void setName(TString name) {fChanName = name;}
  void setNumber(Int_t number) {fChanNumber=number;}
  void setConnectedDevice(TString cd) {fConnectedDevice=cd;}
  void setVoltageMin(Float_t vmin) {fVoltageMin=vmin;}
  void setVoltageMax(Float_t vmax) {fVoltageMax=vmax;}
  void setVoltageRange(Float_t vmin, Float_t vmax) {fVoltageMin=vmin; fVoltageMax=vmax;}
  void setNTriggers(int nt);

  TString getName(){return fChanName; }
  Int_t   getChanNumber(){return fChanNumber; }
  TString getConnectedDevice() {return fConnectedDevice;}
  Float_t getVoltageMin() {return fVoltageMin; }
  Float_t getVoltageMax() {return fVoltageMax; }
  Float_t getVoltageStep() { return -1;}
  int getNTriggers() const;

private:

  void updateHardwareSpecs();

  Int_t     fChanNumber;  //Channel number (0, 1, 2, ...)
  TString   fChanName;    //Channel name (e.g. "CH_A")
  TString   fConnectedDevice; //Name of device connected to this channel (e.g. "PMT_0")
  Float_t   fVoltageMin;  //Min voltage on this channel
  Float_t   fVoltageMax;  //Max voltage on this channel
  int   _nTriggers;


  ClassDef(ScopeDataInfo,1)
};

#endif
