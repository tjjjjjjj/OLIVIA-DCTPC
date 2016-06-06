#ifndef DMTPC_DAQ_HH
#define DMTPC_DAQ_HH
#include "TROOT.h"


class MaxCamCamera;
class MaxCamConfig;
class TH2F;
class TH1F;
class TDatime;
class TTimeStamp;
class TTree;
class TString;
class TFile;
class TCanvas;
class TPad;
class MaxCamTrack;
#include <vector>
using std::vector;
#include <map>
using std::map;
class ScopeHandler;
#include "DmtpcDataset.hh"
#include "ScopeTypes.hh"
//#include "MaxCamCamera.hh"
//#include "ScopeHandler.hh"
#include "MaxCamChannel.hh"


class DmtpcDAQ  {

public:

  // Constructors
  DmtpcDAQ() {}
  DmtpcDAQ(int debugLevel, TString filePath, TString cameraType, TString scopeType);

  virtual ~DmtpcDAQ() {};

  virtual int beginRun();
  virtual int beginRunCV1();
  virtual int endRun();
  
  virtual int makeBiasFrames(int nimages=100,bool zeroExposure=false);

  virtual int beforeEvent();
  virtual int event();
  virtual int afterEvent();

  virtual const char* GetName() { return "DmtpcDAQ"; }

  void saveEvent();
  void clearEventMemory();

  virtual bool isCCDTriggered(int icam);
  virtual bool isScopeTriggered(int isco);
  
  void saveFile();
  void setSaveFlag(bool doSave=false) { _doSave=doSave; }
  bool getSaveFlag() { return _doSave; }
  const char* getFileName();

  void setAutoSave(bool autoSave=false) { _autoSave=autoSave; }
  bool getAutoSave() { return _autoSave; }

  void fillTree(); 

  // camera:
  unsigned int nCCD() { return _camList.size(); }
  vector<MaxCamCamera*>  *camList() { return &_camList; }
  vector<MaxCamCamera*> _camList;
  MaxCamCamera* ccd(int i) { return _camList[i]; }

  // digitizer:
  ScopeHandler* scope() { return _scopeHandler; }
  unsigned int nScope();
  
  //     digitizer trigger
  void setScopeTrigLevel(int iScope, int itrg, unsigned int level);
  void setScopeTrigSlope(int iScope, int itrg, ScopeTriggerSlope slope);
  void setScopeTrigSource(int iScope, int itrg, ScopeTriggerSource source);
  void setScopeTrigEngine(int iScope, int itrg, ScopeTriggerEngine engine);
  void setScopeTrigEngineOperation(int iScope, ScopeTriggerEngOp operation);
  unsigned int getScopeTrigLevel(int iScope, int itrg);
  ScopeTriggerSlope getScopeTrigSlope(int iScope, int itrg);
  ScopeTriggerSource getScopeTrigSource(int iScope, int itrg);
  ScopeTriggerEngine getScopeTrigEngine(int iScope, int itrg);
  ScopeTriggerEngOp getScopeTrigEngineOperation(int iScope);

  float getScopeTriggerLevel(int iScope, int iChannel);
  unsigned int getScopeTriggerLevel1(int iScope);
  unsigned int getScopeTriggerLevel2(int iScope);
  void  setScopeTriggerLevel1(int iScope, unsigned int level);
  void  setScopeTriggerLevel2(int iScope, unsigned int level);

  float getScopeSamplingRate(int iScope);
  void  setScopeSamplingRate(int iScope);

  unsigned int getScopeVoltageRange(int iScope, int iChannel);
  void  setScopeVoltageRange(int iScope, int iChannel, unsigned int range); 
  float getScopeVoltageMin(int iScope, int iChannel);
  float getScopeVoltageMax(int iScope, int iChannel);
  void  setScopeInputImpedance(int iScope, int iChannel, unsigned int impedance); 
  float getScopeInputImpedance(int iScope, int iChannel);

  void  setScopeInputCoupling(int iScope, int iChannel, ScopeInputCoupling coupling); 
  ScopeInputCoupling getScopeInputCoupling(int iScope, int iChannel);

  
  unsigned int  setScopeTriggerDelay(int iScope, unsigned int delay); 

  unsigned int setScopeTriggerSlope( int iScope, unsigned int slope);

  //void setScopeConfigState(int state);
  //int  getScopeConfigState();

  DmtpcDataset* data() { return _data; }

  int setGlobalExposureTime(long msec);

  virtual TString getDBAccessFileName() {return _dbAccessFile;}
  virtual void setDBAccessFileName(const char * filename)
    {if (filename) _dbAccessFile = filename;}

protected:

 MY_MYSQL * getDBHandle() {return db_handle;} 

private:

  ScopeHandler *_scopeHandler;

  DmtpcDataset *_data;

  double _deltaTemp; // deviation from desired CCD temperature

  bool _doSave; // Automatically save event to file (otherwise use saveEvent method)

  bool _autoSave; // save tree such that others can view it in progress 
  
  TString _scopeType;

  MY_MYSQL * db_handle;

  TString _dbAccessFile;

  ClassDef(DmtpcDAQ,0)
};

#endif

