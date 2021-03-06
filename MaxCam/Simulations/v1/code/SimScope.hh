#ifndef __SIMSCOPE__
#define __SIMSCOPE__

//ROOT includes
#include "TMath.h"
#include "TObject.h"
#include "TString.h"

//Viper includes
#include "ScopeDataInfo.hh"
#include "../../../ScopeWaveformData.hh"

class SimScope : public TObject {

public: 

  SimScope();
  SimScope(TString ser, Int_t pmtnum);
  //SimScope(const SimScope &other);
  virtual ~SimScope();

  void setSerialNumber(TString ser){fSerialNumber = ser;}
  void setScopeNumber(Int_t num){fScopeNumber = num;}
  void setBoardType(TString bt);
  void setClockRate(Float_t cr) {fClockRate = cr;}
  void setRecordPreSize(Int_t sz) {fRecordPreSize = sz;}
  void setRecordLength(Int_t sz) {fRecordLength = sz;}
  
  //TJ's functions
  void setDecayTime(double Td){fDecayTime = Td;}
  void setRiseTime(double Tr){fRiseTime = Tr;}
  void setDriftVel(double vd){fDriftVel = vd;}
  void setScale(double Scale){fScale = Scale;}
  void setScopeRes(double scoperes){fScopeRes = scoperes;}
  void setVoltageNoise(double SD){fVoltageNoise = SD;}
  void setTimeNoise1(double a1){fTimeNoise1 = a1;}
  void setTimeNoise2(double a2){fTimeNoise2 = a2;}
  void setZeroOffset(int zero_offset){fZeroOffset = zero_offset;}
  void setZeroOffsetVar(double zero_offset_var){fZeroOffsetVar=zero_offset_var;}
  void setNSamples(int Nsamples){fNSamples = Nsamples;}
  void setVerticalOffset(double vertical_offset){fVerticalOffset = vertical_offset;}
  void setDecayOffset(int decay_offset){fDecayOffset = decay_offset;}
  void setAttenuation(int atten){fAttenuation = atten;}

  double getAttenuation(){return fAttenuation;}
  double getDecayTime(){return fDecayTime;}
  double getRiseTime(){return fRiseTime;}
  double getDriftVel(){return fDriftVel;}
  double getScale(){return fScale;}
  double getScopeRes(){return fScopeRes;}
  double getVoltageNoise(){return fVoltageNoise;}
  double getTimeNoise1(){return fTimeNoise1;}
  double getTimeNoise2(){return fTimeNoise2;}
  int getZeroOffset(){return fZeroOffset;}
  double getZeroOffsetVar(){return fZeroOffsetVar;}
  int getNSamples(){return fNSamples;}
  double getVerticalOffset(){return fVerticalOffset;}
  int getDecayOffset(){return fDecayOffset;}
  //end TJ's functions

  ScopeDataInfo* chan(Int_t channum);
  TObjArray* chans() { return fChans; }
  TString getSerialNumber(){return fSerialNumber;}
  Int_t   getScopeNumber(){return fScopeNumber;}
  TString getBoardType(){return fBoardType;}
  Float_t getClockRate() {return fClockRate;}
  Int_t   getRecordLength() { return fRecordLength; }
  Int_t   getRecordPreSize() { return fRecordPreSize; }
  Int_t   getNbitsPerSamp() { return fNbitsPerSamp;}
  Int_t   getNlevels() {return (Int_t)(TMath::Power(2,getNbitsPerSamp()));}
  Float_t getVoltageStep(Int_t chan);
  Float_t getLevelOfVolts(Int_t chan, Float_t level_volts);

  ScopeWaveformData* wf(Int_t itrig);	
  void setWfs() { setWfs(fScope, fScopeChan);}
  void setWfs(SimScope* scope, Int_t chan);
  void resetWfs();
  void assignScope(SimScope* scope, Int_t chan) {fScope=scope; fScopeChan=chan;}
  
  
private:

  void updateHardwareSpecs();

  TObjArray* fChans;      //Array of SimScopeChan 
  TString fSerialNumber;  //Scope Serial Number
  Int_t   fScopeNumber;   //Scope Number
  TString fBoardType;     //Type of scope board --> defines nchannels, and others
  Int_t   fNChan;         //Number of channels on this scope
  Int_t   fNbitsPerSamp;  //Bits per sample (=8 for Alazar ATS860)
  SimScope*   fScope;        //Pointer to scope that this PMT plugs into
  Int_t       fScopeChan;    //Scope channel number that this PMT plugs into
  Float_t fClockRate;     //Digitization rate in GHz (max is 0.25 for ATS860)
  Int_t   fRecordPreSize; //Number of samples in the pre-trigger waveform
  Int_t   fRecordLength;  //Total number of samples per waveform
  TObjArray*  fWfs;         //Array of waveforms in an event, one per trigger

  double fRiseTime;
  double fDecayTime;
  double fDriftVel;
  double fScale;
  double fScopeRes;
  double fTimeNoise1;
  double fTimeNoise2;
  double fVoltageNoise;
  int fZeroOffset;
  double fZeroOffsetVar;
  int fNSamples;
  double fVerticalOffset;
  int fDecayOffset;
  double fAttenuation;

  ClassDef(SimScope,1)
};

#endif
