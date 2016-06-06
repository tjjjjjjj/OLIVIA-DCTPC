#include "DmtpcDataset.hh"
#include "MaxCamChannel.hh"
//#include "MaxCamImage.hh"
#include "TDatime.h"
#include "TTimeStamp.h"
#include "MaxCamConfig.hh"
#include "TChain.h"
#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TObjString.h"

#ifdef DM_DAQ
#include <mysql/mysql.h>
#endif
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

using namespace std;

ClassImp(DmtpcDataset)


//____________________
//
DmtpcDataset::DmtpcDataset() {
  //

  _event = new DmtpcEvent;
  _comment = new TObjString;
  _location = new TObjString;
  _keyword = new TObjString;
  _detId = new TObjString; 
  _inittime = new TTimeStamp();
  _imageTree = 0;
  _file = 0;
  _listOfDetectorParts = 0;
}

DmtpcDataset::DmtpcDataset(const DmtpcDataset &other) : TObject() {
  cout << GetName() << "Copy Constructor not done" << endl;
  _file = other._file;
  _event = other._event;
  _comment = other._comment;
  _location = other._location;
  _keyword = other._keyword;
  _detId = other._detId; 
  _inittime = other._inittime;
}

DmtpcDataset::~DmtpcDataset() {
  delete _event;
  delete _comment;
  delete _location; 
  delete _keyword;
  delete _detId;
  delete _inittime;

  if (_imageTree) delete _imageTree;
  
  if (_file){
     if (_file->IsOpen()) _file->Close();
     delete _file;
  }
  if (_listOfDetectorParts) delete _listOfDetectorParts;
}


void
DmtpcDataset::fill() { 
  _imageTree->Fill();  
}

void
DmtpcDataset::autoSave() {
  _imageTree->AutoSave("SaveSelf");
}


void
DmtpcDataset::getEvent(int i) {
  // this looks like a typo
  // should be || instead of &&... (jb 10/28/09)
    if (i<0 && i>=tree()->GetEntries()) {
        cout << GetName() << ": event number incorrect! " << i << endl;
        return;
    }
    tree()->GetEvent(i);
}

TH2F*
DmtpcDataset::getBiasFrame(int iccd) {
    TString bfname="biasFrame";
    bfname+=iccd;
    TH2F* bframe=(TH2F*)file()->Get(bfname);
    if (!bframe) {
        cout << GetName() << ": no bias frame for camera " << iccd
             << "  with name: " << bfname
             << endl;
    }
    return bframe;
}

TH2F*
DmtpcDataset::getBiasFrameOverscan(int iccd) {
    TString bfname="biasFrameOverscan";
    bfname+=iccd;
    TH2F* bframe=(TH2F*)file()->Get(bfname);
    if (!bframe) {
        cout << GetName() << ": no bias frame for camera " << iccd
             << "  with name: " << bfname
             << endl;
    }
    return bframe;
}


Int_t 
DmtpcDataset::getNcameras() {
  return getNbiasFrames();
}


Int_t 
DmtpcDataset::getNbiasFrames() {
  
  Int_t iccd = 1;  // biasFrame1 is the lowest-numbered bias frame
  while (1) {
    TString bfname = "biasFrame";
    bfname += iccd;
    TH2F* bframe=(TH2F*)file()->Get(bfname);
    if (!bframe) {
      iccd--;
      break;
    }
    if (iccd > 10) break; // just in case...
    iccd++;
  }

  return iccd;
}

void 
DmtpcDataset::createRootFile( const char *fname, TString foption) {


    foption.ToLower();    

    if (_file){
      if (_file->IsOpen()) _file->Close();
      delete _file;
    }
    if (_imageTree) delete _imageTree;
    _imageTree = 0;
    //
    // write
    //
    if (fname && foption=="recreate") {
      _file = TFile::Open(fname, foption);
      _imageTree = new TTree("dmtpc",  "DMTPC events");
      tree()->Branch("event", "DmtpcEvent", &_event, 32000, 0);
    }
    //
    // read
    //
    else if (fname && foption=="") {
        _file = TFile::Open(fname);
        _imageTree = new TChain("dmtpc",  "DMTPC events");
        chain()->Add(fname);
        chain()->SetBranchAddress("event", &_event);
        if(_file->FindKey("keyword") != NULL)
           _keyword=(TObjString*)_file->Get("keyword");
        if(_file->FindKey("location") != NULL)
           _location=(TObjString*)_file->Get("location");
        if(_file->FindKey("comment") != NULL)
           _comment=(TObjString*)_file->Get("comment");
    }
    else if (fname) assert(0);

}

const char* DmtpcDataset::getFileName() { return _file->GetName(); }


void 
DmtpcDataset::clearEventMemory() {
  // remove data from TClonesArray 
  event()->ccdData()->Delete(); 
  event()->overscan()->Delete();
  event()->ccdConfig()->Delete();
  event()->scopeData()->Delete();
  event()->scopeDataInfo()->Delete();
  event()->experimentConfig()->Delete();
  event()->mcTrack()->Delete();
  event()->mcCcdDigi()->Delete();

//  event()->~DmtpcEvent(); 
}


void 
DmtpcDataset::openRootFile( const char *fname) {
  createRootFile(fname,"");
}

void
DmtpcDataset::closeRootFile(){
  if (_imageTree) delete _imageTree;
  _imageTree = 0;
  if (!_file) return;
  if (!_file->IsOpen()) return;
  _file->Close();
  delete _file; 
  _file = 0; 
}




int
DmtpcDataset::setRunNumberFromDB(const char *dbAccessFile) {
#ifdef DM_DAQ

  char server[256];
  char user[256];
  char pass[256];
  char database[256];
  ifstream infile(dbAccessFile);
  infile >> server >> user >> pass >> database;

  // get name of the last file in the DB
  MYSQL mysql;
  

  mysql_init(&mysql);
  if (!mysql_real_connect( &mysql, server, user, pass, database,0, NULL, 0 )) {
    fprintf(stderr, "Failed to connect to database:  Error: %s\n", mysql_error(&mysql));
    cout << GetName() << ": Cannot connect to DB" << endl;
    return 0;
  }

  TString id = _detId->String(); 

  cout<< GetName() << ": detector id = " << id << endl;

  TString sql("SELECT * FROM run_desc"); 
  if (id!="") {
    sql+=TString(" WHERE detid = '");
    sql+=id;
    sql+=TString("'");
  }

  sql+=TString(" ORDER BY timestamp DESC LIMIT 1");
  mysql_real_query(&mysql, (const char*)sql, sql.Length());

  cout << GetName() << ": sql=" << sql << endl;

  MYSQL_RES *res=mysql_store_result(&mysql);

  if (!res || !mysql_num_rows(res)) {
    cout << GetName() << ": No entries, setting run number to 1" << endl;
    event()->setRunNumber(1);
    return 0;
  }

  MYSQL_ROW  row=mysql_fetch_row(res);
  cout << GetName() << ": last run name " << row[1] << endl;

  int runno=0;
  char runfn[100]; 
  if (id == "") strcpy(runfn,"dmtpc_run"); 
  else sprintf(runfn,"dmtpc_%s_",id.Data());  
  strcat(runfn, "%d.root"); 

  sscanf(row[1],runfn, &runno); // field 1 is run name: 
  mysql_close(&mysql);
  runno++; // choose next available run number

  cout <<GetName()<<": New run number = " << runno << endl;
  event()->setRunNumber(runno);

#else

  cout << GetName() << ": Cannot get run number - MySQL access not active! " << endl
       << dbAccessFile << " not used" << endl;
  
#endif
  return 0;
}



void
DmtpcDataset::setComment(TString c) {
  *_comment=(const char*)c;
}

void
DmtpcDataset::setDetId(TString c)
{
  *_detId = (const char*)c;
}

void
DmtpcDataset::setKeyword(TString c) {
  *_keyword=(const char*)c;
}



void
DmtpcDataset::setLocation(TString c) {
  *_location=(const char*)c;
}


int
DmtpcDataset::saveRun(char *scpdest, char *fname) {
  // Uploads run file to storage directory and fills a database with run information.

#ifdef DM_DAQ

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
    cout << GetName() << ": Cannot connect to DB" << endl;
    return 0;
  }


  // find new name
  int runno=event()->runNumber();
  TString zeros="0000";
  if (runno>9999) zeros="";
  else if (runno>999) zeros="0";
  else if (runno>99)  zeros="00";
  else if (runno>9)   zeros="000";
  TString newName = TString("dmtpc_"); 
  TString id = _detId->String(); 
  if (id=="") newName+=TString("run");
  else 
  {
    newName+=id; 
    newName+=TString("_"); 
  }
  newName += TString(zeros);
  newName += runno;
  newName += TString(".root");


  // find userid
  const char *usid = gSystem->Getenv("USER");

  //find exposure time (asking db is weird, but asking an event's MaxCamConfig doesn't work!)
/*  int exposureTime;
  TString sqlex = "SELECT exposure, ccdid FROM ccd ORDER BY timestamp DESC LIMIT 1";
  mysql_real_query(&mysql, (const char*)sqlex, sqlex.Length());
  MYSQL_RES * res = mysql_store_result(&mysql);
  if(res == NULL || !mysql_num_rows(res))
    {
      exposureTime=0;
    }
  else
    {
      MYSQL_ROW row = mysql_fetch_row(res);
      exposureTime = atoi(row[0]); 
    }

      

  mysql_free_result(res); 


  //calculate average variables
  double anode_i_avg, anode_i_rms, temp0_avg, temp0_rms; 
  double pressure_avg, pressure_rms, anode_hv_avg, anode_hv_rms; 
  double drift_hv_avg, drift_hv_rms, humidity_avg, humidity_rms;

  TString sqlavg = "SELECT avg(value), std(value) FROM wire_i WHERE timestamp>'";
  sqlavg+=_inittime->AsString("s");
  sqlavg+="'";


  cout << sqlavg << endl;
  
  mysql_real_query(&mysql, (const char*)sqlavg, sqlavg.Length());
  res = mysql_store_result(&mysql);
  if(res == NULL || !mysql_num_rows(res))
    {anode_i_avg=0; anode_i_rms=0;}
  else
    {
      MYSQL_ROW row = mysql_fetch_row(res);
      anode_i_avg = atof(row[0]); anode_i_rms=atof(row[1]);
    }

  mysql_free_result(res); 

  sqlavg.ReplaceAll("wire_i","wire_hv");

  cout << sqlavg << endl;


  mysql_real_query(&mysql, (const char*)sqlavg, sqlavg.Length());
  res = mysql_store_result(&mysql);
  if(res == NULL || !mysql_num_rows(res))
    {anode_hv_avg=0; anode_hv_rms=0;}
  else
    {
      MYSQL_ROW row = mysql_fetch_row(res);
      anode_hv_avg = atof(row[0]); anode_hv_rms=atof(row[1]);
    }

  mysql_free_result(res); 

  sqlavg.ReplaceAll("wire_hv","mesh_hv");

  cout << sqlavg << endl;

  mysql_real_query(&mysql, (const char*)sqlavg, sqlavg.Length());
  res = mysql_store_result(&mysql);
  if(res == NULL || !mysql_num_rows(res))
    {drift_hv_avg=0; drift_hv_rms=0;}
  else
    {
      MYSQL_ROW row = mysql_fetch_row(res);
      drift_hv_avg = atof(row[0]); drift_hv_rms=atof(row[1]);
    }

  mysql_free_result(res);

  sqlavg.ReplaceAll("mesh_hv","pressure");

  TString sqlp = "SELECT value from pressure WHERE timestamp>'";
  sqlp +=  _inittime->AsString("s");
  sqlp+="'";

  cout << sqlp << endl;
  mysql_real_query(&mysql, (const char*)sqlp, sqlp.Length());
  MYSQL_RES * resp = mysql_store_result(&mysql);
  if(mysql_num_rows(resp)>0)
    {
      cout << sqlavg << endl;
      
      mysql_real_query(&mysql, (const char*)sqlavg, sqlavg.Length());
      res = mysql_store_result(&mysql);
      if(res == NULL || !mysql_num_rows(res))
	{pressure_avg=0; pressure_rms=0;}
      else
	{
	  MYSQL_ROW row = mysql_fetch_row(res);
	  if(row[0]!="" && row[1]!="")
	    {pressure_avg = atof(row[0]); pressure_rms=atof(row[1]);}
	  else
	    pressure_avg=0; pressure_rms=0;
	}

      mysql_free_result(res);  
    }
  else
    {
      mysql_free_result(resp);
      sqlp = "SELECT value FROM pressure ORDER BY timestamp DESC LIMIT 1";
      mysql_real_query(&mysql, (const char*)sqlp, sqlp.Length());
      resp = mysql_store_result(&mysql);
      MYSQL_ROW row = mysql_fetch_row(resp);
      pressure_avg = atof(row[0]); pressure_rms=0;

      mysql_free_result(resp);
    }
  
  sqlavg.ReplaceAll("pressure","temp0");

    cout << sqlavg << endl;

  mysql_real_query(&mysql, (const char*)sqlavg, sqlavg.Length());
  res = mysql_store_result(&mysql);
  if(res == NULL || !mysql_num_rows(res))
    {temp0_avg=0; temp0_rms=0;}
  else
    {
      MYSQL_ROW row = mysql_fetch_row(res);
      temp0_avg = atof(row[0]); temp0_rms=atof(row[1]);
    }

  mysql_free_result(res);  

  sqlavg.ReplaceAll("temp0","humidity");
  cout << sqlavg << endl;

  mysql_real_query(&mysql, (const char*)sqlavg, sqlavg.Length());
  res = mysql_store_result(&mysql);
  if(res == NULL || !mysql_num_rows(res))
    {humidity_avg=0; humidity_rms=0;}
  else
    {
      MYSQL_ROW row = mysql_fetch_row(res);
      humidity_avg = atof(row[0]); humidity_rms=atof(row[1]);
    }

  mysql_free_result(res);


  //The following is highly specific to the number of cameras, and the addition of saving the serial numbers to the DB

  double ccd_temp_avg_top, ccd_temp_rms_top; 
  double ccd_avgpixel_avg_top, ccd_avgpixel_rms_top;
  double ccd_temp_avg_bottom, ccd_temp_rms_bottom; 
  double ccd_avgpixel_avg_bottom, ccd_avgpixel_rms_bottom;

  sqlavg = "SELECT avg(temperature),std(temperature) FROM ccd WHERE timestamp>'";
  sqlavg += _inittime->AsString("s");
  sqlavg += "' && ccdid='A80333'";

  cout << sqlavg << endl;
  
  mysql_real_query(&mysql, (const char*)sqlavg, sqlavg.Length());
  res = mysql_store_result(&mysql);
  if(res == NULL || !mysql_num_rows(res))
    {ccd_temp_avg_top=0; ccd_temp_rms_top=0;}
  else
    {
      MYSQL_ROW row = mysql_fetch_row(res);
      if(row != NULL)
	{
	  unsigned long* leng = mysql_fetch_lengths(res);
	  if(leng[0] !=0)
	    ccd_temp_avg_top = atof(row[0]); 
	  else
	    ccd_temp_avg_top = 0;
	  if(leng[1]!=0)
	    ccd_temp_rms_top=atof(row[1]);
	  else
	    ccd_temp_rms_top =0;
	}
    }

  mysql_free_result(res);

  sqlavg.ReplaceAll("temperature","avgpixel");

  cout << sqlavg << endl;
  
  mysql_real_query(&mysql, (const char*)sqlavg, sqlavg.Length());
  res = mysql_store_result(&mysql);
  if(res == NULL || !mysql_num_rows(res))
    {ccd_avgpixel_avg_top=0; ccd_avgpixel_rms_top=0;}
  else
    {
      MYSQL_ROW row = mysql_fetch_row(res);
      unsigned long* leng = mysql_fetch_lengths(res);
      if(leng[0]!=0)
	ccd_avgpixel_avg_top = atof(row[0]); 
      else
	ccd_avgpixel_avg_top =0;
      if(leng[1] != 0)
	ccd_avgpixel_rms_top=atof(row[1]);
      else
	ccd_avgpixel_rms_top = 0;
    }

  mysql_free_result(res);

  sqlavg.ReplaceAll("avgpixel","temperature");
  sqlavg.ReplaceAll("A80333","A80334");

  cout << sqlavg << endl;
      
  mysql_real_query(&mysql, (const char*)sqlavg, sqlavg.Length());
  res = mysql_store_result(&mysql);
  if(res == NULL || !mysql_num_rows(res))
    {ccd_temp_avg_bottom=0; ccd_temp_rms_bottom=0;}
  else
    {
      MYSQL_ROW row = mysql_fetch_row(res);
      unsigned long* leng = mysql_fetch_lengths(res);
      if(leng[0] !=0)
	ccd_temp_avg_bottom = atof(row[0]); 
      else
	ccd_temp_avg_bottom = 0;
      if(leng[1] !=0)
	ccd_temp_rms_bottom=atof(row[1]);
      else
	ccd_temp_rms_bottom =0;
    }

  mysql_free_result(res);
  
  sqlavg.ReplaceAll("temperature","avgpixel");
  cout << sqlavg << endl;

  mysql_real_query(&mysql, (const char*)sqlavg, sqlavg.Length());
  res = mysql_store_result(&mysql);
  if(res == NULL || !mysql_num_rows(res))
    {ccd_avgpixel_avg_bottom=0; ccd_avgpixel_rms_bottom=0;}
  else
    {
      MYSQL_ROW row = mysql_fetch_row(res);
      unsigned long* leng = mysql_fetch_lengths(res);
      if(leng[0] !=0)
	ccd_avgpixel_avg_bottom = atof(row[0]);
      else
	ccd_avgpixel_avg_bottom = 0;
      if(leng[1] != 0)
	ccd_avgpixel_rms_bottom=atof(row[1]);
      else
	ccd_avgpixel_rms_bottom = 0;
    }

  mysql_free_result(res);
*/
  // copy file
  TString copyCommand="cp ";
  if (TString(scpdest).Contains(":")) copyCommand="scp ";
  TString scp = copyCommand + getFileName() + " " + scpdest + newName;


  // fill DB
  TString sql="INSERT INTO run_desc (file, keyword, location, description, timestamp, userid, detid";

#ifdef DMTPC_BADBOYS
  sql += ", exposure, anode_i_avg, anode_i_rms, anode_hv_avg, anode_hv_rms, drift_hv_avg, drift_hv_rms, pressure_avg, pressure_rms, temp0_avg, temp0_rms, humidity_avg, humidity_rms, ccd_temp_avg_top, ccd_temp_rms_top, ccd_temp_avg_bottom, ccd_temp_rms_bottom, ccd_avgpixel_avg_top, ccd_avgpixel_rms_top, ccd_avgpixel_avg_bottom, ccd_avgpixel_rms_bottom";
#endif

  sql += ") VALUES( '";
  sql += newName;
  sql += TString("', '");
  sql += _keyword->String();
  sql += TString("', '");
  sql += _location->String();
  sql += TString("', '");
  sql += _comment->String();
  sql += TString("', NOW(), '");
  sql += usid;
  sql += TString("',");
  sql+=TString("'");
  sql += id; 
#ifdef DMTPC_BAD_BOYS
  sql+=TString("', '");
  sql+=exposureTime;
  sql+=TString("', '");
  sql+=anode_i_avg;
  sql+=TString("', '");
  sql+=anode_i_rms;
  sql+=TString("', '");
  sql+=anode_hv_avg;
  sql+=TString("', '");
  sql+=anode_hv_rms;
  sql+=TString("', '");
  sql+=drift_hv_avg;
  sql+=TString("', '");
  sql+=drift_hv_rms;
  sql+=TString("', '");
  sql+=pressure_avg;
  sql+=TString("', '");
  sql+=pressure_rms;
  sql+=TString("', '");
  sql+=temp0_avg;
  sql+=TString("', '");
  sql+=temp0_rms;
  sql+=TString("', '");
  sql+=humidity_avg;
  sql+=TString("', '");
  sql+=humidity_rms;
  sql+=TString("', '");
  sql+=ccd_temp_avg_top;
  sql+=TString("', '");
  sql+=ccd_temp_rms_top;
  sql+=TString("', '");
  sql+=ccd_temp_avg_bottom;
  sql+=TString("', '");
  sql+=ccd_temp_rms_bottom;
  sql+=TString("', '");
  sql+=ccd_avgpixel_avg_top;
  sql+=TString("', '");
  sql+=ccd_avgpixel_rms_top;
  sql+=TString("', '");
  sql+=ccd_avgpixel_avg_bottom;
  sql+=TString("', '");
  sql+=ccd_avgpixel_rms_bottom;
#endif
  sql +=   "')";
		   
  cout << sql << endl;

  mysql_real_query(&mysql, (const char*)sql, sql.Length());
  mysql_close(&mysql);

  cout << GetName() << ": copying file to destination " << endl;
  cout << GetName() << ": " << scp << endl;

  gSystem->Exec(scp);

  _savedLocation=TString(scpdest)+newName;

  cout << GetName() << ": file finished copying to " << _savedLocation << endl;
  //cout << "file finished copying to " << TString(scpdest)+newName << endl;

#else

  cout << GetName() << ": MySQL access not active! " << endl
       << scpdest << "  " << fname << endl;
#endif

  return 0;
}
