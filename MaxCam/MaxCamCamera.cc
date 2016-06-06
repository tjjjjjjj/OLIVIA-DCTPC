#include "MaxCamCamera.hh"
#include "MaxCamConfig.hh"
#include "DmtpcDataConverter.hh"
#include "TH2.h"
#include "TH1.h"
#include "TFile.h"

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
using std::flush;
using namespace DmtpcDataConverter;
#include <fstream>

#include <mysql/mysql.h>



ClassImp(MaxCamCamera)

//____________________
//  Base class for camera daq.
//   
MaxCamCamera::MaxCamCamera()  {
  config = new MaxCamConfig("ccdconfig","CCD Camera configuration");;
  img_histo = new TH2S;
  _sw = new TStopwatch;
  _biasFrame=0;
  _biasFrameOverscan=0;
  _debugFlag=0;
}


void
MaxCamCamera::deleteBiasFrame() {
  // Deletes the bias image so subsequent images are not
  // corrected for the bias.
  if (_biasFrame) delete _biasFrame;
  if(_biasFrameOverscan) delete _biasFrameOverscan;
  _biasFrame=0;
  _biasFrameOverscan=0;
}



int
MaxCamCamera::makeBiasFrame(int nimages, bool zeroExposure) {
  // Take 100 bias frames (no exposure, no shutter) and make one
  // average bias frame. As a result, all subsequent images will have 
  // the bias frame subtracted. In order to stop the bias subtraction
  // call deleteBiasFrame().
  
  cout << endl << GetName() << ": Make bias frame: ";

  deleteBiasFrame();

  setDarkFrame(); 
  long expotime=config->exposureTime;
  if (zeroExposure) setExposureTime(0); 

  TH2F* _biasFrameOverscantmp=0;

  for (int i=0; i<nimages; i++) {
    clickCamera();
    if (!i) {
      TString biasFrameName="biasFrame";
      biasFrameName+=getConfiguration()->cameraID;
      _biasFrame=ccdExpand(createHisto( biasFrameName));
      if(getConfiguration()->digitizeOverscan)
	{
	  TString biasFrameOverscanName="biasFrameOverscan";
	  biasFrameOverscanName+=getConfiguration()->cameraID;
	  _biasFrameOverscantmp = ccdExpand(createOverscanHisto(biasFrameOverscanName));
	}
    } 
    else { 
      TH2S *tmpimage = createHisto("addbias");
      _biasFrame->Add( tmpimage );
      delete tmpimage;
      if(getConfiguration()->digitizeOverscan)
	{
	  TH2S* tmposimage = createOverscanHisto("addosbias");
	  _biasFrameOverscantmp->Add(tmposimage);
	  delete tmposimage;
	}
    }
    cout << "." << flush;
  }
  cout << endl;

  TString biasFrameOverscanName="biasFrameOverscan";
  biasFrameOverscanName+=getConfiguration()->cameraID;
  if(_biasFrameOverscantmp) 
    {
      _biasFrameOverscan=(TH2F*)_biasFrameOverscantmp->Clone(biasFrameOverscanName);
      delete _biasFrameOverscantmp;
    }


  if (_biasFrame) _biasFrame->Scale(1./nimages);
  if (_biasFrameOverscan) _biasFrameOverscan->Scale(1./double(nimages));

  TString bfile="biasFrame_";
  bfile += config->cameraID;
  bfile += ".root";

  TFile fbias(bfile,"RECREATE");
  if (_biasFrame) _biasFrame->Write();  
  if (_biasFrameOverscan) _biasFrameOverscan->Write();
  fbias.Write();
  fbias.Close();

  // restore exposure time if reset to zero for bias frame
  // (nonzero for dark frame)
  setExposureTime(expotime);
  // restore shutter operation
  setNormalFrame();

  return 0;
}



int
MaxCamCamera::readBiasFrame(TString fname) {
  // Read bias frame from a file
  
  deleteBiasFrame();
  clickCamera();
  TString bname="biasFrame";
  bname+=config->cameraID;
  _biasFrame=ccdExpand(createHisto(bname));
  _biasFrame->Reset();
  if(getConfiguration()->digitizeOverscan)
    {
      TString bosname="biasFrameOverscan";
      bosname+=config->cameraID;
      _biasFrameOverscan=ccdExpand(createOverscanHisto(bosname));
      _biasFrameOverscan->Reset();
    }

  if (fname="") {
    fname = "biasFrame_";
    fname += config->cameraID;
    fname += ".root";
  }
  cout << GetName() << ": Reading bias frame " << fname << endl;

  TFile biasfile(fname);
  TString biasFrameName="biasFrame";
  biasFrameName+=config->cameraID;
  TH2F* tmpframe=(TH2F*)biasfile.Get(biasFrameName);
  assert (tmpframe);
  _biasFrame->Add(tmpframe);
  if(getConfiguration()->digitizeOverscan)
    {      
      TString bosname="biasFrameOverscan";
      bosname+=config->cameraID;
      TH2F* tmposframe = (TH2F*)biasfile.Get(bosname);
      assert(tmposframe);
      _biasFrameOverscan->Add(tmposframe);
    }
  
  biasfile.Close();

  return 0;
}




int
MaxCamCamera::setExposureTime(long msec) {
  config->exposureTime=msec;
  return 0;
}

int
MaxCamCamera::setDigitizeOverscan(bool digi) { 
    cout << GetName() << ": digitization of overscan not implemented " << digi << endl;
    return -1;
}


void
MaxCamCamera::writeCCDConfigToDB(const char *fname, unsigned int * db_handle) {

  MYSQL * mysql; 
  MYSQL  mysql_struct; 

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
      cout << "MaxCamChannel: Cannot connect to DB" << endl;
      return;
    }
  }
  else
  {
    mysql = (MYSQL*)db_handle;  
    mysql_ping(mysql); 
  }

  // fill DB
  TString sql=TString("INSERT INTO ccd (temperature, exposure, daqtime, avgpixel, ccdid) VALUES( '");

  sql += config->CCDTemp;
  sql += TString("', '");
  sql += config->exposureTime;
  sql += TString("', '");
  sql += config->daqTime;
  sql += TString("', '");
  sql += float(config->meanForExposedPixels);
  sql += TString("', '");
  sql += config->serialNumber;
  sql += TString("' )");

  mysql_real_query(mysql, (const char*)sql, sql.Length());
  if (!db_handle)
    mysql_close(mysql);
  
}
