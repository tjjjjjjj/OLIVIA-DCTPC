#ifndef DMTPC_Run_HH
#define DMTPC_Run_HH

#include "DmtpcDAQ.hh"

class TMapFile;
class DmtpcEvent;

class DmtpcRun : public DmtpcDAQ  {

public:

  DmtpcRun(int debugLevel=0, TString outputFile="test.root", TString ccdType="dummy", TString scopeType="dummy");

  virtual ~DmtpcRun() {};

  virtual int beforeEvent();
  virtual int afterEvent();

  virtual const char* GetName() { return "DmtpcRun"; }

  static void setDriftVoltage(double val, const char *table="drift_hv");

  static void setAnodeVoltage(double val, const char *table="anode_hv");

  static void setPressureAndSave(double val);

  static void rampUp();

  static void rampDown();


private:

  //TMapFile *_DqmFile;
  //DmtpcEvent *_DqmEvent;

  void saveCCDInfoToDB();
  void saveScopeInfoToDB();
  
  ClassDef(DmtpcRun,0)
};

#endif

