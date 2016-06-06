#include "MaxCamTest.hh"
#include "MaxCamConfig.hh"
//#include "MaxCam.hh"
#include "MaxCamAlta.hh"
#include "MaxCamChannel.hh"
#include "MaxCamTrack.hh"

#include <iostream>
#include <fstream>
#include <vector>
using std::cout;
using std::endl;
using std::cerr;
using std::vector;
#include <mysql/mysql.h>

#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TDatime.h"
#include "TFile.h"
#include "TPad.h"
#include "TTree.h"
#include "TSystem.h"
#include "TNtuple.h"

#include<map>
using std::map;
using std::pair;

ClassImp(MaxCamTest)

//____________________
// Class to run MaxCam camera. 
//
//
  MaxCamTest::MaxCamTest(int debugLevel, TString outputFile, TString ctype) : MaxCamTestBase(debugLevel, outputFile, ctype) {
  // Constructor initializes a camera (one for now) and
  // creates a tree for data-taking.

  _eventCounter= new int[2];
  _eventCounter[0]=0;
  _eventCounter[1]=0;
    
}




int
MaxCamTest::event() {
  
  for (unsigned int i=0; i<_camList.size(); i++) ccd(i)->getTemperature();
  dtime()->Set();

  wire()->readFromDB("wire_hv");

  mesh()->readFromDB("mesh_hv");

  pressure()->readFromDB("pressure");

  int nt=0;
  bool isOk=false;
  while(1) {
    _eventCounter[0]++;

    // this really works up to 2 cameras
    // otherwise, make separate processes
    for (unsigned int i=0; i<_camList.size(); i++) {
      ccd(i)->expose();
    }
    for (unsigned int i=0; i<_camList.size(); i++) {
      ccd(i)->grabImage ();
    }

    writeCCDConfigToDB();

    if ( ++nt>getTriggerTrials() ) break;

    for (unsigned int i=0; i<_camList.size(); i++) {
      if (isTriggered(i)) { isOk=true; break; }
    }
    if (isOk) break;
  }

  if (getSaveFlag()) {
    saveEvent();
  }

  if (isOk) _eventCounter[1]++;

  return 0;
}


int
MaxCamTest::end() {
  // Ends datataking and closes root file.

  // save counter
  TNtuple *nt=new TNtuple("counter","","tot:acc");
  nt->Fill( _eventCounter[0], _eventCounter[1] );
  nt->Write();

  MaxCamTestBase::end();

  cout << "            BUT, IF YOU WANT TO PERMANENTLY KEEP THE RUN, PLEASE USE saveRun !" << endl;

  return 0;
}




int
MaxCamTest::saveRun(char *keyword, char *description, char *location, char *scpdest, char *fname) {
  // Uploads run file to storage directory and fills a database with run information.


  char server[256];
  char user[256];
  char pass[256];
  char database[256];
  ifstream infile(fname);
  infile >> server >> user >> pass >> database;

  // get name of the last file in the DB
  MYSQL mysql;
  
  mysql_init(&mysql);
  if (!mysql_real_connect( &mysql, server, user, pass, database,0, NULL, 0 )) {
    cout << "MaxCamChannel: Cannot connect to DB" << endl;
    return 0;
  }

  TString sql("SELECT * FROM run_desc ORDER BY timestamp DESC LIMIT 1");
  mysql_real_query(&mysql, (const char*)sql, sql.Length());

  //cout << sql << endl;

  MYSQL_RES *res=mysql_store_result(&mysql);

  if (!mysql_num_rows(res)) {
    cout << "MaxCamChannel: No entries" << endl;
    return 0;
  }

  MYSQL_ROW  row=mysql_fetch_row(res);

  int runno=0;
  sscanf(row[1],"ccdrun_%d.root", &runno); // field 1 is run name: 

  // find new name
  runno++;
  TString zeros="0000";
  if (runno>9999) zeros="";
  else if (runno>999) zeros="0";
  else if (runno>99)  zeros="00";
  else if (runno>9)   zeros="000";
  TString newName = TString("ccdrun_") + zeros;
  newName += runno;
  newName += TString(".root");

  // copy file
  TString scp = TString("scp ") + getFileName() + " " + scpdest + newName;

  // fill DB
  sql=TString("INSERT INTO run_desc (file, keyword, location, description, timestamp) VALUES( '");

  sql += newName;
  sql += TString("', '");
  sql += keyword;
  sql += TString("', '");
  sql += location;
  sql += TString("', '");
  sql += description;
  sql += TString("', NOW() )");
		   
  cout << sql << endl;

  mysql_real_query(&mysql, (const char*)sql, sql.Length());
  mysql_close(&mysql);
  
  cout << "copying file to destination" << endl;

  gSystem->Exec(scp);

  return 0;
}


void 
MaxCamTest::writeValueToDB(TString table, double val) {
  int iread=0, iwrite=0;
  if (table=="wire_hv")       { iread=1; iwrite=1; }
  else if (table=="pressure") { iread=2; iwrite=-1; }
  MaxCamChannel ch((const char*)table,"", iread,iwrite,"dbaccess.txt");
  ch.readFromDB((const char*)table);
  ch.setValue = val;
  ch.saveValueToDBAndCheck((const char*)table);
}


void
MaxCamTest::writeCCDConfigToDB(const char *fname) {

  char server[256];
  char user[256];
  char pass[256];
  char database[256];
  ifstream infile(fname);
  infile >> server >> user >> pass >> database;

  MYSQL mysql;
  
  mysql_init(&mysql);
  if (!mysql_real_connect( &mysql, server, user, pass, database,0, NULL, 0 )) {
    cout << "MaxCamChannel: Cannot connect to DB" << endl;
    return;
  }

  TH1F *hy=drawYields(0);
  float avgPixel = hy->GetMean()/(ccd()->config->row_width*ccd()->config->img_rows);

  // fill DB
  TString sql=TString("INSERT INTO ccd (temperature, exposure, daqtime, avgpixel) VALUES( '");

  sql += ccd()->config->CCDTemp;
  sql += TString("', '");
  sql += ccd()->config->exposureTime;
  sql += TString("', '");
  sql += ccd()->config->daqTime;
  sql += TString("', '");
  sql += avgPixel;
  sql += TString("' )");
		   
  mysql_real_query(&mysql, (const char*)sql, sql.Length());
  mysql_close(&mysql);
  
}


void
MaxCamTest::recoverFromDischarge(UInt_t sleepTime) {
  MaxCamChannel mesh("mesh_hv","", 0,0,"dbaccess.txt");
  mesh.readFromDB("mesh_hv");
  double mesh0=mesh.setValue;

  MaxCamChannel wire("wire_hv","", 1,1,"dbaccess.txt");
  wire.readFromDB("wire_hv");
  double wire0=wire.setValue;

  // ramp down
  cout << "MaxCamTest: RAMP DOWN and sleep for " << sleepTime << "sec" << endl;
  writeValueToDB("mesh_hv", 0.02);
  writeValueToDB("wire_hv", 0.02);

  gSystem->Sleep(sleepTime); // 15*60 sec 

  
  // ramp up
  cout << "MaxCamTest: RAMP UP" << endl;
  writeValueToDB("mesh_hv", mesh0);
  writeValueToDB("wire_hv", wire0);

  gSystem->Sleep(25000); // 
}


void
MaxCamTest::scanVoltage(float min, float max, float step, bool wire) {

  setTriggerTrials(0);
  setSaveFlag(true);

  for (int i=0; i<10; i++) {

    for (float v=min; v<max; v+=step) {
      cout << "Setting wire voltage to " << v << "kV" << endl;
      for (int dummy=0; dummy<10; dummy++) {
        gSystem->Sleep(200);
        if (wire) writeValueToDB("wire_hv",v); 
        else writeValueToDB("mesh_hv",v);
      }
      gSystem->Sleep(1000);
      for (int j=0; j<10; j++) { event(); drawImage(0, "colz"); gPad->Update(); }
    }

  }


  

}
