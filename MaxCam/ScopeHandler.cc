// ScopeHandler.cc
//
#include <iostream>
//#include <string>
#include <vector>
#include <assert.h>

#include "ScopeHandler.hh"
#include "ScopeAlazarATS860.hh"
#include "ScopeData.hh"
#include "ScopeDataInfo.hh"

using std::cout;
using std::endl;
//using std::string;
using std::vector;
#include <fstream>

#include <mysql/mysql.h>

ScopeHandler::~ScopeHandler() {}

ScopeHandler::ScopeHandler() {}

//ScopeHandler::ScopeHandler(const string* scopeType) {
ScopeHandler::ScopeHandler(const TString* scopeType) {

  // determine type of scope
  // Instantiate scope object
  // Add scope to list of scopes
  // Open communication w/ scope
  // Do some default setup (maybe, maybe not)

  U32 nscope=AlazarBoardsFound();
  cout << GetName() << ": constructor" << endl;
  cout << GetName() << ": scopeType = " << *scopeType << endl;
  cout << GetName() << ": found boards = " << nscope << endl;

  for (U32 isc=0; isc< nscope; isc++)  addScope(scopeType);
}

float
ScopeHandler::getMaxVoltage(int ii, U32 voltRangeSpecifier) { return scope(ii)->getMaxVoltage(voltRangeSpecifier); }

//float
//ScopeHandler::getDtOfClockRate(int ii, AU32 sampleRate) { return scope(ii)->getDtOfClockRate(sampleRate); }

int ScopeHandler::getNChannels(int isco) { return scope(isco)->getNChannels(); }
int ScopeHandler::getNLevels(int isco) { return scope(isco)->getNLevels(); }

int
ScopeHandler::acquireTriggers(int ii, float duration_ms) {
  // acquire triggers from scope ii for duration_ms (milliseconds)

  return scope(ii)->acquireTriggers(data(ii), duration_ms);
}

int
ScopeHandler::readData(int ii) {
  // acquire triggers from scope ii for duration_ms (milliseconds)

  return scope(ii)->readData(data(ii));
}

int
ScopeHandler::getNValidTriggers(int ii) { return data(ii)->getNValidTriggers(); }

vector<ScopeWaveformData*>
ScopeHandler::getWaveforms(int ii) { return data(ii)->getWaveforms(); }

// this might need to be PRIVATE...
ScopeConfig*
ScopeHandler::getScopeConfig(int ii) {
  return scope(ii)->getScopeConfig();
}

float
ScopeHandler::getSamplingRate(int isco) { return scope(isco)->getSamplingRate(); }

float
ScopeHandler::getVoltageMin(int isco, int chanNum) { return scope(isco)->getVoltageMin(chanNum); }
float
ScopeHandler::getVoltageMax(int isco, int chanNum) { return scope(isco)->getVoltageMax(chanNum); }

ScopeInputCoupling
ScopeHandler::getInputCoupling(int isco, int chanNum) { return scope(isco)->getInputCoupling(chanNum); }

float
ScopeHandler::getInputImpedance(int isco, int chanNum) { return scope(isco)->getInputImpedance(chanNum); }

float
ScopeHandler::getVoltageStep(int isco, int chanNum) { return scope(isco)->getVoltageStep(chanNum); }

float
ScopeHandler::getTriggerLevel(int isco, int chanNum) { return scope(isco)->getTriggerLevel(chanNum); }

ScopeTriggerSlope
ScopeHandler::getTriggerSlope(int isco, int chanNum) { return scope(isco)->getTriggerSlope(chanNum); }


void 
ScopeHandler::addScope(const TString* scopeType) {
//ScopeHandler::addScope(const string* scopeType) {

  validateScopeType(scopeType);

  if (*scopeType == Scope::ALAZAR_ATS860) {
    int boardID=_scopeList.size();
    cout << "adding scope: " << *scopeType << " to list, with boardID=" << boardID << endl;
    _scopeList.push_back( new ScopeAlazarATS860(boardID) );
  } else {
    cout << scopeType << " is not yet implemented" << endl;
  }

  // keep a list of the scope types in use
  setScopeType(scopeType);
}

void
ScopeHandler::openScope(int ii) {
  scope(ii)->openScope();
}

void
ScopeHandler::openScope() {
  U32 nscope=AlazarBoardsFound();
  for (U32 isc=0; isc< nscope; isc++)  scope(isc)->openScope();
}


ScopeConfig*
ScopeHandler::configureDAQState(int ii, Scope::DAQ_State daqstate) {
  cout << GetName() << ".configureDAQState() " << endl;
  return scope(ii)->configureDAQState(daqstate);
}


//////////////////////////////////////////////////////
// set/get trigger methods using Alazar-specific types...
void
ScopeHandler::setTrigLevel(int isco, int itrg, unsigned int level) { scope(isco)->setTrigLevel(itrg, level); }
unsigned int 
ScopeHandler::getTrigLevel(int isco, int itrg) { return scope(isco)->getTrigLevel(itrg); }

void 
ScopeHandler::setTrigSlope(int isco, int itrg, ScopeTriggerSlope slope) { scope(isco)->setTrigSlope(itrg, slope); }
ScopeTriggerSlope
ScopeHandler::getTrigSlope(int isco, int itrg) { return scope(isco)->getTrigSlope(itrg); }

void 
ScopeHandler::setTrigSource(int isco, int itrg, ScopeTriggerSource source) { scope(isco)->setTrigSource(itrg, source); }
ScopeTriggerSource
ScopeHandler::getTrigSource(int isco, int itrg) { return scope(isco)->getTrigSource(itrg); }

void 
ScopeHandler::setTrigEngine(int isco, int itrg, ScopeTriggerEngine engine) { scope(isco)->setTrigEngine(itrg, engine); }
ScopeTriggerEngine
ScopeHandler::getTrigEngine(int isco, int itrg) { return scope(isco)->getTrigEngine(itrg); }

void 
ScopeHandler::setTrigEngineOperation(int isco, ScopeTriggerEngOp operation) { scope(isco)->setTrigEngineOperation(operation); }
ScopeTriggerEngOp 
ScopeHandler::getTrigEngineOperation(int isco) { return scope(isco)->getTrigEngineOperation(); }


// end set/get trigger 
//////////////////////////////////////////////////////





void
ScopeHandler::setTriggerLevel1(int ii, U32 level) {
  return scope(ii)->setTriggerLevel1(level);
}
U32
ScopeHandler::getTriggerLevel1(int ii) {
  return scope(ii)->getTriggerLevel1();
}
void
ScopeHandler::setTriggerLevel2(int ii, U32 level) {
  return scope(ii)->setTriggerLevel2(level);
}
U32
ScopeHandler::getTriggerLevel2(int ii) {
  return scope(ii)->getTriggerLevel2();
}

void
ScopeHandler::setTriggerLevel1Volts(int ii, int chanNumber, float level_volts) {
  return scope(ii)->setTriggerLevel1Volts(chanNumber, level_volts);
}
void 
ScopeHandler::setVoltageRange(int scopeNumber, int chanNumber, U32 range) { scope(scopeNumber)->setVoltageRange(chanNumber, range); }

U32
ScopeHandler::getVoltageRange(int scopeNumber, int chanNumber) { return scope(scopeNumber)->getVoltageRange(chanNumber); }

void 
ScopeHandler::setInputImpedance(int scopeNumber, int chanNumber, U32 impedance) { scope(scopeNumber)->setInputImpedance(chanNumber, impedance); }
void 
ScopeHandler::setInputCoupling(int scopeNumber, int chanNumber, ScopeInputCoupling coupling) { scope(scopeNumber)->setInputCoupling(chanNumber, coupling); }

U32
ScopeHandler::setTriggerDelay(int scopeNumber, unsigned int delay) { return scope(scopeNumber)->setTriggerDelay(delay); }

U32
ScopeHandler::setTriggerSlope(int scopeNumber, unsigned int slope) { return scope(scopeNumber)->setTriggerSlope(slope); }


void 
ScopeHandler::addScopeData(int ii, ScopeConfig* sc) {
  _scopeDatasets.push_back( new ScopeData(sc) );
}

void
ScopeHandler::addScopeDataInfo(int ii) { _scopeDataInfo.push_back( new ScopeDataInfo() ); }

ScopeData*
ScopeHandler::data(int ii) { return _scopeDatasets.at(ii); }

ScopeDataInfo*
ScopeHandler::dataInfo(int ii) { return _scopeDataInfo.at(ii); }

void
ScopeHandler::clearData(int ii) { data(ii)->clearData(); }

//
// PRIVATE METHODS
//

Scope* 
ScopeHandler::scope(int ii) { return _scopeList[ii]; }

void 
//ScopeHandler::validateScopeType(const string* scopetype) {
ScopeHandler::validateScopeType(const TString* scopetype) {
  if (!Scope::isValidType(scopetype)) {
    assert(!"invalid scope type");
  }
}

//const string* 
const TString* 
ScopeHandler::getScopeType(int ii) { return _scopeTypes.at(ii); }

void
//ScopeHandler::setScopeType(const string* scopetype) {
ScopeHandler::setScopeType(const TString* scopetype) {
  _scopeTypes.push_back(scopetype);
}

char*
ScopeHandler::GetName() { return "ScopeHandler"; }





void
ScopeHandler::saveScopeInfoToDB(const char *fname, MYSQL * db_handle) {

  int npulse[10][2];
  float epulse[10][2];
  for (unsigned int isco=0; isco<nScope(); isco++) {
    for (int ich=0; ich<getNChannels(isco); ich++) {
      npulse[isco][ich]=0;
      epulse[isco][ich]=0;
    }
  }
  
  for (unsigned int isco=0; isco<nScope(); isco++) {
    for (int ich=0; ich<getNChannels(isco); ich++) {

      unsigned int nwfch=data(isco)->dataChan(ich)->wfs()->size();
      for (unsigned int i=0; i<nwfch; i++) {
	TH1 *hwf=data(isco)->dataChan(ich)->wf(i);
	npulse[isco][ich]++;
	epulse[isco][ich]+=hwf->GetMaximum();
      }

    }
  }


  MYSQL * mysql; 
  MYSQL mysql_struct; 

  if (db_handle == NULL)
  {
    char server[256];
    char user[256];
    char pass[256];
    char database[256];
    ifstream infile(fname);
    infile >> server >> user >> pass >> database;
  
    mysql = &mysql_struct;
    mysql_init(mysql);
    if (!mysql_real_connect( mysql, server, user, pass, database,0, NULL, 0 )) {
      cout << GetName() <<": Cannot connect to DB" << endl;
      return;
    }
  }
  else
  { 
    mysql = db_handle; 
    mysql_ping(mysql); 
  }
  

      
  for (unsigned int isco=0; isco<nScope(); isco++) {
    for (int ich=0; ich<getNChannels(isco); ich++) {

      // fill DB
      TString sql=TString("INSERT INTO scope (scopeid, chid, nwf, esum) VALUES( '");
      
      sql += isco;
      sql += TString("', '");
      sql += ich;
      sql += TString("', '");
      sql += npulse[isco][ich];
      sql += TString("', '");
      sql += epulse[isco][ich];
      sql += TString("' )");
      
      //cout << sql << endl;

      mysql_real_query(mysql, (const char*)sql, sql.Length());

      if (!db_handle)
        mysql_close(mysql);
      
    }

  }
  
}
