
#include "DmtpcRun.hh"
#include "MaxCamChannel.hh"
#include "MaxCamCamera.hh"
#include "MaxCamConfig.hh"
#include "ScopeHandler.hh"

//#include "TMapFile.h"
#include "TH2F.h"

#include <iostream>
#include <fstream>
#include <vector>
using std::cout;
using std::endl;
using std::cerr;
using std::vector;



ClassImp(DmtpcRun)


//____________________
// Class to collect data for the DMTPC experiment. 
// 
// It's based on DmtpcDAQ, but includes database support, data 
// quality and archiving.
//
//
//
DmtpcRun::DmtpcRun(int debugLevel, TString outputFile, TString ccdType, TString scopeType) :
    DmtpcDAQ(debugLevel, outputFile, ccdType, scopeType) {

    //_DqmFile =  TMapFile::Create("dmtpc_dqm.root", "RECREATE", 100000, "DQM file for DMTPC run.");
    //_DqmEvent = new DmtpcEvent;

    data()->file()->cd();
}


int
DmtpcRun::beforeEvent() {
  // Prepare event.

  //
  // read values that are less likely to change from event-to-event
  //

  //TStopwatch sw;
  //sw.Start();


  // temperature from CCD camera
  for (unsigned int i=0; i<_camList.size(); i++) ccd(i)->getTemperature();
  
  // voltages, currents and pressure from database
  int nPars=0;
  
    MaxCamChannel* driftHV   = new MaxCamChannel("driftHV",    "Drift Voltage", 0,  0, "dbaccess.txt");
  driftHV->readFromDB("drift_hv",false,getDBHandle());
  new( (*data()->event()->experimentConfig())[nPars++] ) MaxCamChannel(*driftHV);
  delete driftHV;
  
  MaxCamChannel* anodeHV   = new MaxCamChannel("anodeHV",    "Anode Voltage", 1,  1, "dbaccess.txt");
  anodeHV->readFromDB("anode_hv",false,getDBHandle());
  new( (*data()->event()->experimentConfig())[nPars++] ) MaxCamChannel(*anodeHV);
  delete anodeHV;

  MaxCamChannel* driftI    = new MaxCamChannel("driftI",     "Drift current", 2, -1, "dbaccess.txt");
  driftI->readFromDB("drift_i",false,getDBHandle());
  new( (*data()->event()->experimentConfig())[nPars++] ) MaxCamChannel(*driftI);
  delete driftI;

  MaxCamChannel* temp0   = new MaxCamChannel("temp0",    "Temperature point 0", 4,  -1, "dbaccess.txt");
  temp0->readFromDB("temp0",false,getDBHandle());
  new( (*data()->event()->experimentConfig())[nPars++] ) MaxCamChannel(*temp0);
  delete temp0;

  MaxCamChannel* anodeI    = new MaxCamChannel("anode_i",     "Anode current", 5, -1, "dbaccess.txt");
  if (anodeI) { 
    anodeI->readFromDB("anode_i",false,getDBHandle());
    new( (*data()->event()->experimentConfig())[nPars++] ) MaxCamChannel(*anodeI);
    delete anodeI;
  }

  MaxCamChannel* pressure = new MaxCamChannel("pressure_cdg",   "Gas pressure",  -1, -1, "dbaccess.txt");    
  if (pressure) {
    pressure->readFromDB("pressure",false,getDBHandle());
    new( (*data()->event()->experimentConfig())[nPars++] ) MaxCamChannel(*pressure);
    delete pressure;
  }

  //
  // set values that change every event
  //
  DmtpcDAQ::beforeEvent();

  //cout << "Before Event Time: " << sw.RealTime() << endl;

  return 0;
}



int
DmtpcRun::afterEvent() {
  // After event.
  
  // dispatch event for DQM
  //*_DqmEvent=*data()->event();
  //_DqmFile->Update();

  //TStopwatch sw;
  //sw.Start();

  cout << GetName() << "->afterEvent() " << endl;
  saveCCDInfoToDB();
  cout << GetName() << "->afterEvent() post CCDInfoToDB() " << endl;
  saveScopeInfoToDB();
  cout << GetName() << "->afterEvent() post ScopeInfoToDB() " << endl;

  // clear memory
  int ret= DmtpcDAQ::afterEvent();
  cout << GetName() << "->afterEvent() post DmtpcDAQ::afterEvent() " << endl;

  //cout << "After Event Time: " << sw.RealTime() << endl;

  return ret;
}


void 
DmtpcRun::setDriftVoltage(double val, const char *table) {
  MaxCamChannel mesh(table,"", 0,0,"dbaccess.txt");
  mesh.readFromDB(table);
  mesh.setValue = val;
  mesh.saveValueToDBAndCheck(table);
}

void 
DmtpcRun::setAnodeVoltage(double val, const char *table) {
  MaxCamChannel wire(table,"", 1,1,"dbaccess.txt");
  wire.readFromDB(table);
  wire.setValue = val;
  wire.saveValueToDBAndCheck(table);
}


void 
DmtpcRun::setPressureAndSave(double val) {
  MaxCamChannel pressure("pressure","", 3,-1,"dbaccess.txt");
  pressure.readFromDB("pressure");
  pressure.setValue = val;
  pressure.saveValueToDBAndCheck("pressure");
}

void 
DmtpcRun::rampUp() {
  MaxCamChannel hvstatus("hvstatus","", -1,-1,"dbaccess.txt");
  hvstatus.setValue = 1;
  hvstatus.saveValueToDBAndCheck("hvstatus");
}

void 
DmtpcRun::rampDown() {
  MaxCamChannel hvstatus("hvstatus","", -1,-1,"dbaccess.txt");
  hvstatus.setValue = -1;
  hvstatus.saveValueToDBAndCheck("hvstatus");
}


void
DmtpcRun::saveCCDInfoToDB() {
  for (unsigned int i=0; i<_camList.size(); i++) {
    ccd(i)->writeCCDConfigToDB("dbaccess.txt", (unsigned int *) getDBHandle());
  }
}


void
DmtpcRun::saveScopeInfoToDB() {
  // Scope data is saved into database for on-line data-quality monitoring. 
  scope()->saveScopeInfoToDB("dbaccess.txt", (MYSQL*)getDBHandle());
}
