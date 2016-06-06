#include "DmtpcDAQ.hh"
#include "MaxCamConfig.hh"

#ifdef CAM_FLI
#include "MaxCam.hh"
#endif

#ifdef CAM_ALTA
#include "MaxCamAlta.hh"
#endif

#ifdef CAM_ANDOR
#include "MaxCamAndor.hh"
#include "atmcdLXd.h"
#endif

#include "MaxCamDummyCamera.hh"
#include "MaxCamCamera.hh"
#include "MaxCamChannel.hh"
#include "MaxCamTrack.hh"
#include "MaxCamImageTools.hh"
#include "Scope.hh"
#include "ScopeData.hh"
#include "ScopeDataInfo.hh"
#include "ScopeWaveform.hh"
#include "ScopeHandler.hh"
//#include "ScopeConfig.hh"
#include "ScopeTypes.hh"
#include "MaxCamImage.hh"
#include "DmtpcEvent.hh"
#include "DmtpcDataConverter.hh"

using namespace DmtpcDataConverter;

#include <iostream>
#include <fstream>
#include <vector>
using std::cout;
using std::endl;
using std::cerr;
using std::vector;

#include "TObjString.h"
#include "TString.h"
#include "TFile.h"
#include "TH2.h"
#include "TDatime.h"
#include "TTimeStamp.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStopwatch.h"
#include "TRandom.h"

//#include "TThread.h"
//#include <stdlib.h>

#include<map>
using std::map;
using std::pair;


ClassImp(DmtpcDAQ)

//____________________
// Class to collect data for the DMTPC experiment. 
// 
// Normally, one would use the following sequence:
//
// beginRun(); // search and initialize hardware, prepare data output
// 
// // collect events - repeat this for every event
// beforeEvent();
// event();
// afterEvent();
//
// endRun(); // close data file and hardware handles
//
//
  DmtpcDAQ::DmtpcDAQ(int debugLevel, TString outputFile, TString cameraType, TString scopeType)  {
  // Constructor initializes CCD cameras and oscilloscopes, and
  // creates a data structure for output.

  // Set the default db access filename here
  setDBAccessFileName("dbaccess.txt");

  //
  //   Camera initialization 
  //
  cameraType.ToLower();
  if (cameraType=="dummy") {

    for (int i=0; i<2; i++) { // add 2 fake cameras
      cout << GetName() << ":  DUMMY CCD CAMERA" << endl;
      _camList.push_back( new MaxCamDummyCamera( debugLevel) );

    }
  }
#ifdef CAM_ALTA
  else if (cameraType=="apogee") {

    MaxCamAlta::discoverCameras();
  
    for (unsigned short i=0; i<MaxCamAlta::nCamera; i++) 
      _camList.push_back( new MaxCamAlta() );
  }
#endif

#ifdef CAM_FLI
  else if (cameraType=="fli") {

    _camList.push_back( new MaxCam(debugLevel) );
  }
#endif

#ifdef CAM_ANDOR
  else if (cameraType=="andor") {

    long int lNumCameras;
    GetAvailableCameras(&lNumCameras);
    for (unsigned short i=0; i<lNumCameras; i++) 
      _camList.push_back( new MaxCamAndor() );
  }
#endif

  else assert(0);
  
  for (unsigned int i=0; i<_camList.size(); i++) {
    assert (!_camList[i]->openCamera(i+1)); 
    // some default setup
    ccd(i)->setExposureTime(0); // exposure
  }
  //_cam->setNumberOfFlushes(16); // cleaning of CCD before exposure
  //_cam->startBkgFlush(); // start flushing
  _deltaTemp=2;

  // 
  //    Digitizer initialization 
  //
  _scopeHandler = new ScopeHandler(&Scope::ALAZAR_ATS860);
  // scope() returns _scopeHandler ...
  scope()->openScope();
  _scopeType = scopeType;

  // 
  //    Event setup
  // 
  _doSave=false;

  //
  //    ROOT output/input
  //
  _data = new DmtpcDataset;
  data()->createRootFile( outputFile, "recreate" );

  db_handle = NULL; 
}

int
DmtpcDAQ::beginRun() {
  // Calibaration of cameras and digitizers.
  cout << GetName() << ": begin run" << endl;

  for (unsigned int i=0; i<_camList.size(); i++) {
    // wait for ccd to cool
    while (ccd(i)->getTemperature()>ccd(i)->getGoalTemperature()+_deltaTemp) {
      cout << GetName() << ": Waiting for CCD " << i << " to cool down, T="<<ccd(i)->getTemperature()
	   <<"  ->  T(goal)=" <<  ccd(i)->getGoalTemperature() << endl;
      gSystem->Sleep(5000);
    }
    // make one fake readout
    ccd(i)->cleanCCD();
  }

  U32 nscope=AlazarBoardsFound();
  for (U32 isc=0; isc < nscope; isc++) {
    scope()->configureDAQState(isc, Scope::DAQ_CHARGE);
    scope()->addScopeData(isc, scope()->getScopeConfig(isc));
    scope()->addScopeDataInfo(isc);
  }
  
  data()->setRunNumberFromDB(getDBAccessFileName()); // use next available number


  char server[256];
  char user[256];
  char pass[256];
  char database[256];
  ifstream infile(getDBAccessFileName());
  infile >> server >> user >> pass >> database;

  db_handle = (MY_MYSQL*) malloc(sizeof(MYSQL)); 
	
  mysql_init((MYSQL*)db_handle); 
  if (!mysql_real_connect((MYSQL*)db_handle, server, user, pass, database, 0, NULL,0))
  {
    fprintf(stderr, "Failed to connect to database:  Error: %s\n", mysql_error((MYSQL*)db_handle));
    cout << GetName() << ": Cannot connect to DB" << endl;
  }

  cout << GetName() << ": Persistent DB Connection created. " << endl ;
  return 0;
}

int DmtpcDAQ::makeBiasFrames(int nimages, bool zeroExposure)
{
  long expotime = ccd(0)->getConfiguration()->exposureTime;
  for(unsigned int i=0; i<_camList.size(); i++)
    {
      ccd(i)->deleteBiasFrame();
      ccd(i)->setDarkFrame();
      if(zeroExposure) ccd(i)->setExposureTime(0);
    }

  vector< TH2F*> _biasFrameOverscantmp;
  vector< TH2F*> _biasFrametmp;

  for(int i=0; i<nimages; i++)
    {
      for(unsigned int j=0; j<_camList.size(); j++)
	ccd(j)->expose();
      
      for(unsigned int j=0; j<_camList.size(); j++)
	{
	  ccd(j)->grabImage();
	  if(!i)
	    {
	      TString biasFrameName = "biasFrame";
	      biasFrameName+=ccd(j)->getConfiguration()->cameraID;
	      TH2S* firstbias = ccd(j)->createHisto(biasFrameName+"S");
	      _biasFrametmp.push_back(ccdExpand(firstbias,biasFrameName));
	      //delete firstbias;
	      firstbias->Delete();
	      if(ccd(j)->getConfiguration()->digitizeOverscan)
		{
		  TString biasFrameOverscanName="biasFrameOverscan";
		  biasFrameOverscanName+=ccd(j)->getConfiguration()->cameraID;
		  TH2S* firstosbias = ccd(j)->createOverscanHisto(biasFrameOverscanName+"S");
		  _biasFrameOverscantmp.push_back(ccdExpand(firstosbias,biasFrameOverscanName));
		  //delete firstosbias;
		  firstosbias->Delete();
		}
	    }
	  else
	    {
	      
	      TH2S* tmpimg = ccd(j)->createHisto("addbias");
	      _biasFrametmp[j]->Add(tmpimg);
	      delete tmpimg;
	      if(ccd(j)->getConfiguration()->digitizeOverscan)
		{
		  TH2S* tmposimage = ccd(j)->createOverscanHisto("addosbias");
		  _biasFrameOverscantmp[j]->Add(tmposimage);
		  delete tmposimage;
		}
	    }
	}
      
      cout << i << "." << flush;
    }
  
  for(unsigned int i=0; i<_camList.size(); i++)
    {
      _biasFrametmp[i]->Scale(1/double(nimages));
      ccd(i)->setBiasFrame(_biasFrametmp[i]);
      if(ccd(i)->getConfiguration()->digitizeOverscan)
	{
	  _biasFrameOverscantmp[i]->Scale(1/double(nimages));
	  ccd(i)->setBiasFrameOverscan(_biasFrameOverscantmp[i]);
	}
      
      TString bfile="biasFrame_";
      bfile += ccd(i)->config->cameraID;
      bfile += ".root";
      
      TFile fbias(bfile,"RECREATE");
      _biasFrametmp[i]->Write();  
      _biasFrameOverscantmp[i]->Write();
      fbias.Write();
      fbias.Close();
      
      // restore exposure time if reset to zero for bias frame
      // (nonzero for dark frame)
      ccd(i)->setExposureTime(expotime);
      // restore shutter operation
      ccd(i)->setNormalFrame();
    }

  return 0;
}

int 
DmtpcDAQ::beginRunCV1() {
  // Calibaration of cameras and digitizers.
  cout << GetName() << ".beginRunCV1()" << endl;
  /*

  for (unsigned int i=0; i<_camList.size(); i++) {
    // wait for ccd to cool
    while (ccd(i)->getTemperature()>ccd(i)->getGoalTemperature()+_deltaTemp) {
      cout << GetName() << ": Waiting for CCD " << i << " to cool down, T="<<ccd(i)->getTemperature()
	   <<"  ->  T(goal)=" <<  ccd(i)->getGoalTemperature() << endl;
      gSystem->Sleep(5000);
    }
    // make one fake readout
    ccd(i)->cleanCCD();
  }

  scope()->configureDAQState(0, Scope::DAQ_CV1);

  scope()->addScopeData(0, scope()->getScopeConfig(0));
  scope()->addScopeDataInfo(0);
  cout << "done with addScopeData" << endl;

  //data()->setRunNumberFromDB(getDBAccessFileName()); // use next available number

  int runNumber = 1;
  data()->event()->setRunNumber(runNumber);
  */
  return 0;
}



bool
DmtpcDAQ::isCCDTriggered(int icam) {
  // Check if a specific CCD needs readout.
  // This needs maping from CCD to scope in the future.
  // Currently, any scope trigger will also trigger
  // all CCD cameras.
  return true;
  //if ( isScopeTriggered(0) )  return true;

  return false;
}


bool
DmtpcDAQ::isScopeTriggered(int isco) {
  // Check if specific oscilloscopes needs readout.
  // Currently, if waveforms are present, event
  // is triggered.

  if ( scope()->data(isco)->getNValidTriggers() > 0 ) return true;

  return false;
}


void DmtpcDAQ::setScopeTrigLevel(int iScope, int itrg, unsigned int level)     { scope()->setTrigLevel(iScope, itrg, level); }
void DmtpcDAQ::setScopeTrigSlope(int iScope, int itrg, ScopeTriggerSlope slope)     { scope()->setTrigSlope(iScope, itrg, slope); }
void DmtpcDAQ::setScopeTrigSource(int iScope, int itrg, ScopeTriggerSource source)   { scope()->setTrigSource(iScope, itrg, source); }
void DmtpcDAQ::setScopeTrigEngine(int iScope, int itrg, ScopeTriggerEngine engine)   { scope()->setTrigEngine(iScope, itrg, engine); }
void DmtpcDAQ::setScopeTrigEngineOperation(int iScope, ScopeTriggerEngOp operation) { scope()->setTrigEngineOperation(iScope, operation); }
unsigned int DmtpcDAQ::getScopeTrigLevel(int iScope, int itrg)  { return scope()->getTrigLevel(iScope, itrg); }
ScopeTriggerSlope DmtpcDAQ::getScopeTrigSlope(int iScope, int itrg)  { return scope()->getTrigSlope(iScope, itrg); }
ScopeTriggerSource DmtpcDAQ::getScopeTrigSource(int iScope, int itrg) { return scope()->getTrigSource(iScope, itrg); }
ScopeTriggerEngine DmtpcDAQ::getScopeTrigEngine(int iScope, int itrg) { return scope()->getTrigEngine(iScope, itrg); }
ScopeTriggerEngOp DmtpcDAQ::getScopeTrigEngineOperation(int iScope)  { return scope()->getTrigEngineOperation(iScope); }





float
DmtpcDAQ::getScopeTriggerLevel(int iScope, int iChannel) {
  // Get trigger level for specified scope. Returns a 
  // floating point value.

  return scope()->getTriggerLevel(iScope, iChannel);
}

unsigned int
DmtpcDAQ::getScopeTriggerLevel1(int iScope) {
  // Get trigger level for specified scope. Returns a 
  // value in ADC units.

  return scope()->getTriggerLevel1(iScope);
}
unsigned int
DmtpcDAQ::getScopeTriggerLevel2(int iScope) {
  // Get trigger level for specified scope. Returns a 
  // value in ADC units.

  return scope()->getTriggerLevel2(iScope);
}
void
DmtpcDAQ::setScopeTriggerLevel1(int iScope, unsigned int level) {
  // Set trigger level for specified scope in ADC units.

  scope()->setTriggerLevel1(iScope, level);
}

void
DmtpcDAQ::setScopeTriggerLevel2(int iScope, unsigned int level) {
  // Set trigger level for specified scope in ADC units.

  scope()->setTriggerLevel2(iScope, level);
}

float 
DmtpcDAQ::getScopeSamplingRate(int iScope) {
  // Get sampling rate

  return scope()->getSamplingRate(iScope);
}

void 
DmtpcDAQ::setScopeSamplingRate(int iScope) {}


unsigned int 
DmtpcDAQ::getScopeVoltageRange(int iScope, int iChannel) {
  // Returns range in ADC units.

  return scope()->getVoltageRange(iScope, iChannel);
}

void 
DmtpcDAQ::setScopeVoltageRange(int iScope, int iChannel, unsigned int range) {
  // Sets scope range in ADC units.

  scope()->setVoltageRange(iScope, iChannel, range);
}

float 
DmtpcDAQ::getScopeVoltageMax(int iScope, int iChannel) {
  // Returns maximum voltage for scope, channel in volts.

  return scope()->getVoltageMax(iScope, iChannel);
}

float
DmtpcDAQ::getScopeVoltageMin(int iScope, int iChannel) {
  // Returns minimum voltage for scope, channel in volts.

  return scope()->getVoltageMin(iScope, iChannel);
}


void 
DmtpcDAQ::setScopeInputImpedance(int iScope, int iChannel, unsigned int impedance) {
  // Sets the input impedance of a channel of the scope

  scope()->setInputImpedance(iScope, iChannel, impedance);
}
float
DmtpcDAQ::getScopeInputImpedance(int iScope, int iChannel) {
  // Returns the input impedance of a channel of the scope
  return scope()->getInputImpedance(iScope, iChannel);
}

void 
DmtpcDAQ::setScopeInputCoupling(int iScope, int iChannel, ScopeInputCoupling coupling) {
  // Sets the input coupling of a channel of the scope
  scope()->setInputCoupling(iScope, iChannel, coupling);
}
ScopeInputCoupling
DmtpcDAQ::getScopeInputCoupling(int iScope, int iChannel) {
  // Returns the input coupling of a channel of the scope

  return scope()->getInputCoupling(iScope, iChannel);
}

unsigned int
DmtpcDAQ::nScope() { return scope()->nScope(); }

unsigned int
DmtpcDAQ::setScopeTriggerDelay( int iScope, unsigned int delay) {

    return scope()->setTriggerDelay(iScope, delay);
}

unsigned int
DmtpcDAQ::setScopeTriggerSlope( int iScope, unsigned int slope) {

    return scope()->setTriggerSlope(iScope, slope);
}


int
DmtpcDAQ::beforeEvent() {
  // Prepare event. Set only experimental parameters
  // that change every event.

  // set time-stamp (precision in seconds)
  data()->event()->timeStamp()->Set();
  data()->event()->UTCtimeStamp()->Set();
  data()->event()->increaseEventNumber();
  // reset save-event-flag
  setSaveFlag(false);

  return 0;
}


////Thread method
//
//void * trigger_cam(void * ccd_ptr)
//{
//  ((MaxCamCamera*)ccd_ptr)->grabImage(); 
//  return (0); 
//}


int
DmtpcDAQ::event() {
  //
  // Calls exposure, trigger collection for every event.
  // Exposures, triggers acquisitions are requested in series.
  // Although start-acquisition loops is fast, small time 
  // misalignment is still possible.
  // 
  // Readout is also done in series, but it takes much longer
  // so this scheme really works up to 2 cameras
  // otherwise it will take too long to download data.
  // In future, plan to make separate processes for camera control.
  //


  //stopwatch for debugging
  //TStopwatch sw;
  //sw.Start();  


  //====================================
  // start acquisition
  //
  unsigned int icam=0;
  // duration in milliseconds (needs to come from scope config.)
  // add time for scanning exposures during ccd readout
  float duration_ms = ccd(0)->config->exposureTime; 
  // -----------> from here
  for (; icam<nCCD(); icam++) ccd(icam)->expose();
  if (_scopeType=="dummy") gSystem->Sleep( UInt_t(duration_ms) );
  else scope()->acquireTriggers(0, duration_ms); // in master-slave operation, only master can trigger (scope==0)
  // ------------ to here, the code should be fast to assure time alignment

  //cout << "Finished acquisition" << endl;

    
  //====================================
  // collect waveforms
  //

  //cout << "Beginning scope readout" << endl;

  unsigned int iisum=0;
  for (unsigned int isco=0; isco<nScope(); isco++) {

    scope()->readData(isco);

    // check for triggeres in each scope
    if (isco==0) { // check if triggered by master board only
      bool isTriggered = isScopeTriggered(isco);
      cout << GetName() << ": " << (isTriggered ? "Master board triggered" : "No trigger") << endl;
      if (!isTriggered) break;
    }

    ////////////////////////////////////
    // Write the number of triggers to the ScopeDataInfo class
    scope()->dataInfo(isco)->setNTriggers(scope()->data(isco)->getNValidTriggers());
    //////////////////////////////////////////////////////////////
    // The Scope configuration will certainly be the same
    // for all triggers coming from the same channel in one event
    // and so we should probably only save the ScopeDataInfo 
    // contents once per channel per event, but since ScopeDataInfo
    // is so small, we do a ScopeDataInfo entry for EVERY trigger.  


    // The TH1F is stored as voltage vs. time 
    // by the acquireTriggers() routine, 

    // pack all waveforms into the TClonesArray
    unsigned int MAX_N_WF = _scopeType!="dummy" ? DmtpcEvent::MAX_N_TRIGGERS : 0;

    for (int ich=0; ich<scope()->getNChannels(isco); ich++) {

    
      //for (unsigned int ii=0; (ii<waveforms.size()) && (ii < MAX_N_WF); ii++ ) {
      unsigned int nwfch=scope()->data(isco)->dataChan(ich)->wfs()->size();
      cout << GetName() << ": Waveforms scope " << isco << " channel " << ich << "............"<<nwfch<<endl;
      for (unsigned int ii=0; (ii<nwfch) && (ii < MAX_N_WF); ii++ ) {
     
	// put a new TH1F into the TClonesArray
	//TString iname = "charge"+ii; // charge0, charge1, etc.
	TString iname = "charge";
	iname += "_";
	iname += isco;
	iname += "_";
	iname += ich;
	iname += "_";
	iname += ii; // charge0, charge1, etc.

	//cout << iname << endl;

	// get the data from the waveform in memory
	ScopeWaveformData* wfChan = scope()->data(isco)->dataChan(ich)->wf(ii);
	//cout << "address of scope waveform = " << wfChan << "  mean = " << wfChan->GetMean() << endl;

	//// put time into underflow bin
	if (!ich) wfChan->SetBinContent( 0, scope()->data(isco)->dataChan(ich)->getTime(ii) ); // trigger time stamp in cha only

	/////////////////////////////////////////////////////////////////
	// put ScopeWaveform into TClonesArray via copy constructor
        // This works since wfChanA is a ScopeWaveform which "is-a" TH1F.
        // Probably want to change this so that the TClonesArray contains
        // ScopeWaveform not TH1F, then you can have ScopeWaveform.timestamp
        // and other auxiliary information about the trace.
        new( (*data()->event()->rawScopeData())[iisum] ) ScopeWaveformData(*wfChan);

	///////////////////////////////////////////////////
	// Save the scope setup into another TClonesArray 
	// in parallel to the waveforms
      	scope()->dataInfo(isco)->setSamplingRate(scope()->getSamplingRate(isco));
      	scope()->dataInfo(isco)->setVoltageRange(scope()->getVoltageMin(isco, ich), scope()->getVoltageMax(isco, ich));
      	scope()->dataInfo(isco)->setChannelId(ich);
      	scope()->dataInfo(isco)->setInputCoupling(scope()->getInputCoupling(isco, ich));
      	scope()->dataInfo(isco)->setInputImpedance(scope()->getInputImpedance(isco, ich));
      	scope()->dataInfo(isco)->setScopeNumber(isco);
      	scope()->dataInfo(isco)->setVoltageStep(scope()->getVoltageStep(isco, ich));
      	scope()->dataInfo(isco)->setTriggerLevel(scope()->getTriggerLevel(isco, ich));
      	scope()->dataInfo(isco)->setTriggerSlope(scope()->getTriggerSlope(isco, ich));
      	//scope()->dataInfo(isco)->setTriggerCause(scope()->getTriggerCause(isco));
	new( (*data()->event()->scopeDataInfo())[iisum] ) ScopeDataInfo(*scope()->dataInfo(isco));

	iisum++;
      }
    }
  } 
  //cout << "Ending scope readout total waveforms saved = " << iisum << endl;
    

  //cout << "Beginning CCD readout" << endl;
  
  //====================================
  // download CCD images
  //
    int iSavedImages=0;

  //Create a bitfield to keep track of which
  //ccd's are triggered
  //char * triggered_ccds = (char*) calloc(nCCD()/8+1,1); 
   
  //vector<TThread*> threads; 

  /* Start the readouts */
  for (unsigned int icam=0; icam<nCCD(); icam++) {
    // check trigger for each camera
    if (!isCCDTriggered(icam)) {
      cout << "ccd not triggered XXXXXXXXXXXXX cancelExposure() " << endl;
      ccd(icam)->cancelExposure();
      continue; 
    }
    cout << GetName() <<": Camera triggered " << icam <<  endl;
    //triggered_ccds[icam/8] |= (1<<(icam%8));
    //while (!ccd(icam)->imageReady());
    setSaveFlag(true); // mark event for saving
    //TThread * th = new TThread(trigger_cam,ccd(icam)); 
    //threads.push_back(th); 
    //}

  /* Loop over threads waiting for them to finish*/
  //for (unsigned int ith = 0; ith < threads.size(); ith++)
  //{
  //  threads[ith]->Join(); 
  //}

  /* Copy data to hists */
  //for (unsigned int icam = 0; icam  > nCCD(); icam++)
    //{
    //if (triggered_ccds[icam/8] & (1 << (icam%8)) == 0) continue; 

    //cout << "Before grabImage()" << endl;
    ccd(icam)->grabImage(); // download image from camera
    

    TString iname="ccd_";
    iname += icam;
    TString oname="overscan_";
    oname+= icam;
    //cout << GetName() <<": Save camera image " << icam <<  endl;
    TH2S *himage=ccd(icam)->createHisto( iname );
    //cout << "Placement new" << endl;
    new( (*data()->event()->rawCcdData())[iSavedImages] ) TH2S(*himage);
    cout << ccd(icam)->getConfiguration()->digitizeOverscan << endl;
    
    if(ccd(icam)->getConfiguration()->digitizeOverscan)
      {
	//cout << "Reached overscan readout" << endl;
	TH2S *osimage=ccd(icam)->createOverscanHisto( oname);
	new( (*data()->event()->rawOverscan())[iSavedImages] ) TH2S(*osimage);
	//cout << "Overscan histo saved" << endl;
	delete osimage;
	//cout << "Overscan histo deleted" << endl;
      }

    new( (*data()->event()->ccdConfig())[iSavedImages++] ) MaxCamConfig(*ccd(icam)->getConfiguration());

    delete himage;
  }
  //free(triggered_ccds); 
  ///* Loop over threads to clear memory */
  //for (unsigned int ith = 0; ith < threads.size(); ith++)
  //{
  //  threads[ith]->Delete(); 
  //}
  //cout << "CCD download time = " << sw.RealTime() << endl;

  //cout << "Ending ccd readout" << endl;



  //====================================
  // save data
  //
  if (getSaveFlag()) {
    saveEvent();
    //data()->tree()->Print();
  }

      
  return 0;
}





int
DmtpcDAQ::afterEvent() {
  // After event.

#ifdef SCOPE_DEBUG

  for (unsigned int isco=0; isco<nScope(); isco++) {
    cout << "afterEvent(): " << scope()->dataInfo(isco)->getSamplingRate() << endl;
  }

  // get info from the TClonesArray (only works if there was a trigger)
  int isco = 0;
  cout << "isScopeTriggered(0) = " << isScopeTriggered(isco) << endl;
  if (isScopeTriggered(isco)) {
    cout << "data()->event()->scopeDataInfo(isco) = "
	 << data()->event()->scopeDataInfo(isco) << endl;

    cout << "data()->event()->scopeDataInfo(isco)->getSamplingRate() = "
	 <<  data()->event()->scopeDataInfo(isco)->getSamplingRate() << endl;
    cout << "data()->event()->scopeDataInfo(isco)->getVoltageMin() = "
	 <<  data()->event()->scopeDataInfo(isco)->getVoltageMin() << endl;
    cout << "data()->event()->scopeDataInfo(isco)->getVoltageMax() = "
	 <<  data()->event()->scopeDataInfo(isco)->getVoltageMax() << endl;
    cout << "data()->event()->scopeDataInfo(isco)->getInputCoupling() = "
	 <<  data()->event()->scopeDataInfo(isco)->getInputCoupling() << endl;
    cout << "data()->event()->scopeDataInfo(isco)->getInputImpedance() = "
	 <<  data()->event()->scopeDataInfo(isco)->getInputImpedance() << endl;
    cout << "data()->event()->scopeDataInfo(isco)->getScopeNumber() = "
	 <<  data()->event()->scopeDataInfo(isco)->getScopeNumber() << endl;
    cout << "data()->event()->scopeDataInfo(isco)->getVoltageStep() = "
	 <<  data()->event()->scopeDataInfo(isco)->getVoltageStep() << endl;
    cout << "data()->event()->scopeDataInfo(isco)->getTriggerLevel() = "
	 <<  data()->event()->scopeDataInfo(isco)->getTriggerLevel() << endl;
    cout << "data()->event()->scopeDataInfo(isco)->getTriggerSlope() = "
	 <<  data()->event()->scopeDataInfo(isco)->getTriggerSlope() << endl;
  }
#endif


  // need to manage the ScopeWaveform data separately:
  // reset scope data buffer.
  for (unsigned int isco=0; isco<nScope(); isco++) 
    scope()->clearData(isco);

  // clear data from event output
  clearEventMemory();

  return 0;
}




void 
DmtpcDAQ::fillTree() {
    // Fill data tree for current event.
    _data->fill();

    //if (_autoSave) _data->autoSave();
}

void
DmtpcDAQ::saveEvent() {
  // Force saving current event into ROOT file.
  fillTree();
}

void 
DmtpcDAQ::clearEventMemory() {
  _data->clearEventMemory();
}


int
DmtpcDAQ::endRun() {
  // Ends datataking and closes root file.

  data()->file()->cd();
  for (unsigned int i=0; i<_camList.size(); i++) {
    if (!ccd(i)->biasFrame()) {
      cout << GetName() << ": NO bias frame for camera " << i << endl;
      continue;
    }
    ccd(i)->biasFrame()->Write();
    if(ccd(i)->getConfiguration()->digitizeOverscan)
      {
	if (!ccd(i)->biasFrameOverscan()) {
	  cout << GetName() << ": NO overscan bias frame for camera " 
	       << i << endl;
	  continue;
	}
      
	ccd(i)->biasFrameOverscan()->Write();
      }
  }

  data()->comment()->Write();
  data()->keyword()->Write();
  data()->location()->Write();
  //data()->listOfDetectorParts()->Write();
 
  saveFile();

  for (unsigned int i=0; i<_camList.size(); i++) ccd(i)->closeCamera();
  cout << GetName() << ": this run has finished and events are saved into " << getFileName() << endl;

  mysql_close((MYSQL*)db_handle); 
  free(db_handle); 
  db_handle = NULL; 

  return 0;
}

void
DmtpcDAQ::saveFile() { data()->write(); data()->file()->Close(); }

const char* 
DmtpcDAQ::getFileName() { return _data->getFileName(); }

int
DmtpcDAQ::setGlobalExposureTime(long msec) {
  // Set exposure time for all CCD cameras and oscilloscopes
  // in the system.

  // ... for CCD's
  for (unsigned int i=0; i<nCCD(); i++) {
    ccd(i)->setExposureTime(msec);       
  }
  // ... for scopes
  for (unsigned int i=0; i<nScope(); i++) {
    cout << GetName() << ": scope acquisition time not set! Using CCD time!"<<endl;
  }

  return 0;
}

