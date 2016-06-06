///////////
// To do //
///////////
// * log data file is not written to after each fill
//   so a viewer would not be able to see recent data...
//   Figure out what triggers a write and do that after each fill
// * currently no timeout for data reads (serial comm, telnet instances, etc)
//   this can block and lead to data loss (I've seen it happen)
// * run serial reads as separate threads?
// * verify automatic log-file roll-over at midnight GMT
// * log-file safety.  if code crashes/restarts, ensure log file is preserved
// * embed the DmtpcLoggableParam within DmtpcLogEntry so that DmtpcLoggableParam doesn't
//   have to be handled/managed by this code
// * communication with database so that DAQ can access certain vars
// * get power status from C++ code, not from external python script...
// * add reality checks on parameter reads. (especially the pressure gauge readings).
// * proper time axes and appropriate labels on TGraphs
// * "AutoSave" option to allow log file to be viewed as it is written
// * need logfile viewing code, both for interactive mode and for automated plot generation
// * could have a single TTimeStamp instance, no?
//   (if the entire program is single-threaded)
//
//
///////////////////////////////////
// Brief description of the code //
///////////////////////////////////
//
//  Code to log the environmental parameters of the 4shooter detector
//
//  A new log file is automatically created every 24 hours at a
//  specified time (currently midnight UTC)
// 
//  Logged data is stored in /data/4shooter/log
//  with filenames of the form:  logYYYYMMDD.root
//
//  The data are all stored in a ROOT tree.  Each tree branch
//  corresponds to a different parameter.  For examle, there is a
//  branch for CF4 pressure, one for the scroll pump power status,
//  etc.  Each entry in a branch is a DmtpcLogEntry instance.
//  DmtpcLogEntry contains a TClonesArray of DmtpcLoggableParam
//  instances.  DmtpcLoggableParam contains a single timestamp and a
//  single value (taken at that time).
//
//  So any given tree entry could contain zero or more readings
//  (DmtpcLoggableParam instances) on each branch.  Each parameter can
//  be read at independent rates (cadences).  The tree filling
//  interval is also independent of the data logging cadences.
//
//  For information on how to stop the logging program, see the
//  comments in the DetectorLogger() function below.
//
// Author:  James Battat  jbattat@mit.edu
// 2010 January 12

//  What needs to be logged (differential vs. single-ended)
//  0  ADC Temperature0
//     ADC Temperature1
//     ADC Temperature2
//     ADC Temperature3
//     ADC Humidity
//  5  ADC TiltmeterX
//     ADC TiltmeterY
//     ADC TiltmeterZ
//     ADC +HV (Anode)
//     ADC -HV (Cathode)
// 10  ADC PMT_HV_0
//     ADC PMT_HV_1
//     ADC PMT_HV_2
//     ADC Big gate valve position indicator
//     ADC Small gate valve position indicator
// 15  ADC 
//     ADC 
//     ADC 
//     ADC 
//     ADC 
// 20  ADC 
//     ADC 
//     ADC 
//     ADC 
//     ADC 
// 25  ADC 
//     ADC 
//     ADC 
//     ADC 
//     ADC 
// 30  ADC 
// 31  ADC 
//
//  
//  CCD Parameters --> but how do you get them?
//    CCD0_temperature   
//    CCD1_temperature
//    CCD2_temperature
//    CCD3_temperature
//
//  Serial Varian Turbo controller state (rpm, temperature)
//  Serial Inficon (CDG, BPG)
//  Serial Convectron
//  Serial MFC (is not logged, but is used)
//
//  telnet power1
//  telnet power2
//  telnet power3
//  telnet power4
//  telnet power5
//  telnet power6
//  telnet power7
//  telnet power8  Scroll pump
//  telnet power9  Big gate
//  telnet power10 Small gate
//  telnet power11 Rough valve
//  telnet power12 MFC valve
//  telnet power13
//  telnet power14 
//  telnet power15
//  telnet power16 
//
//
//


// C/C++ includes
#include <iostream>
#include <sys/time.h>  // gettimeofday()
#include <unistd.h>    // usleep()
#include <fstream>     // ifstream

// ROOT includes
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TDatime.h"
#include "TCanvas.h"

// DMTPC includes
#include "../MaxCamChannel.hh"
#include "../DmtpcLogEntry.hh"
#include "../DmtpcLoggableParam.hh"

using std::cout;
using std::endl;
//using namespace std;

//
// perhaps these #defines could be moved out to an ASCII configuration file
//
#define DMTPC_TPS 1000 // dmtpctime ticks per second
#define SPLIT 0  // splitlevel of the TClonesArray branches
#define MAX_ENTRIES 100  // maximum number of log values between fills
#define USLEEP_US 100000  // sleep time in microseconds

// The Inficon controller (VGC402) reads out two inficon pressure gauges
// The VGC402 has 2 channels.  
// Channel 1 is the BPG (dual gauge, atmosphere to 1e-11 torr)
// Channel 2 is the CDG (capacitance diaphragm gauge, ~1 to 200 torr)
#define INFICON_SERIAL_PORT "/dev/ttyS0"
#define INFICON_BPG_CHAN 1
#define INFICON_CDG_CHAN 2 

Int_t verbose = 1;

typedef Long64_t dmtpctime_t;  // times can be negative

dmtpctime_t logfile_due, // when to start a new log file
  treefill_due = 0,
  temp_due = 0,    // temperatures
  powers_due = 0,  // read power status of IP power strip
  pres_due = 0;    // pressures


dmtpctime_t  dmtpctbase, dmtpcctime;  // start & current epoch times, ms
dmtpctime_t  dmtpctime=0;  // ms since execution started

TFile *logfile;  // the active log file
TTree *logtree;  // the active log file

TString logfiledir("/data/4shooter/log/");
TString logfilename;

/*----------------------------------------------------------------------------
Define logging cadences
----------------------------------------------------------------------------*/
dmtpctime_t treefill_interval = 5*DMTPC_TPS,
  temp_interval  = 2*DMTPC_TPS, 
  power_interval = 3*DMTPC_TPS, 
  pres_interval  = DMTPC_TPS;

/*----------------------------------------------------------------------------
Define variables to hold logging parameters
----------------------------------------------------------------------------*/
// Pressures
TTimeStamp *presTime = new TTimeStamp();  // used by all pressures
MaxCamChannel inficonVGC402_chan("inficonVGC402", "inficonVGC402");  // used by presCF4 and presChamb

DmtpcLogEntry *presCF4 = new DmtpcLogEntry("presCF4", "CDG gauge");
DmtpcLoggableParam *presCF4_param = new DmtpcLoggableParam("presCF4_param", "cdg gauge reading");

DmtpcLogEntry *presChamb = new DmtpcLogEntry("presChamb", "BPG gauge");
DmtpcLoggableParam *presChamb_param = new DmtpcLoggableParam("presChamb_param", "bpg gauge reading");

// Powers
TTimeStamp *powersTime = new TTimeStamp();
DmtpcLogEntry *powerScroll = new DmtpcLogEntry("powerScroll", "Scroll pump power");
DmtpcLoggableParam *powerScroll_param = new DmtpcLoggableParam("powerScroll_param", "Scroll pump power");

//// Temperatures
//TClonesArray *temp_tca_ptr = new TClonesArray("DmtpcLoggableParam", MAX_ENTRIES);
//DmtpcLoggableParam tempParam("temperature0", "cf4 temperature");

// define functions
void do_system_logging();
void do_waiting();
bool stop_requested();
void init_time();
void init_loggable_params();
void init_pressure_params();
dmtpctime_t set_dmtpctime();
dmtpctime_t set_logfiledue();

// data acquisition
void get_pressures();
void get_temperatures();
void get_powers();
void read_powers();
Float_t valOfPowerStatus(TString status);
void fill_tree();
void reset();
void reset_pressures();
void reset_powers();
void reset_temperatures();

/*---------------------------------------------------------------------------
set_dmtpctime():  updates dmtpctime & dmtpcctime for current execution time

Static In/Out dmtpcctime current epoch time, in DMTPC_TPS units
Out dmtpctime  --  time since execution started

All scheduled processing uses 'dmtpctime' to decide when to process.
---------------------------------------------------------------------------*/
dmtpctime_t set_dmtpctime() {

  struct timeval curtime;

  gettimeofday(&curtime, NULL);
  dmtpcctime = (ULong64_t)curtime.tv_sec*DMTPC_TPS + 
    (ULong64_t)curtime.tv_usec*DMTPC_TPS/1000000; // in TICKS
  dmtpctime = dmtpcctime - dmtpctbase; // ticks since execution started

  return dmtpctime;
}

/*---------------------------------------------------------------------------
set_logfiledue(): calc the time from now to the next log file rollover
Static Out  logfile_due  time until log file rolls over (in ticks)
  ---------------------------------------------------------------------------*/
dmtpctime_t set_logfiledue() {
  time_t cursec = time(NULL), logday; // seconds since epoch
  logday = cursec;  // some time within log start day
  // roll at 00:00 next day (rollover is at midnight)
  cout << "logday = " << logday << endl;
  cout << "logday%86400 = " << logday%86400 << endl;
  cout << "dmtpctime = " << dmtpctime << endl;
  logfile_due = (dmtpctime_t)(86400 - logday%86400)*DMTPC_TPS + dmtpctime;
  cout << "logfile_due = " << logfile_due << endl;
  cout << "logfile_due in hours = " << logfile_due/(DMTPC_TPS*3600.) << endl;
  //gmtime_r(&logday, &btime);

  return logfile_due;
}

/*---------------------------------------------------------------------------
begin data logging routines
---------------------------------------------------------------------------*/
void get_temperatures() {
  cout << "   reading temperatures" << endl;
}

void get_powers() {
  cout << "    reading powers" << endl;
  powersTime->Set();
  powerScroll_param->SetTimeStamp(powersTime);
  read_powers();

  powerScroll->Append(powerScroll_param);
}

void get_pressures() {
  cout << "   reading pressures" << endl;

  presTime->Set();  // update the pressure timestamp

  // CDG reading
  presCF4_param->SetValue(inficonVGC402_chan.readInficonController(INFICON_SERIAL_PORT, INFICON_CDG_CHAN));
  presCF4_param->SetTimeStamp(presTime);
  presCF4->Append(presCF4_param);

  // BPG reading
  presChamb_param->SetValue(inficonVGC402_chan.readInficonController(INFICON_SERIAL_PORT, INFICON_BPG_CHAN));
  presChamb_param->SetTimeStamp(presTime);
  presChamb->Append(presChamb_param);

}

void fill_tree() {
  cout << " ------------ fill_tree() -------------- " << endl;
  Int_t nbytes = 0;
  nbytes = logtree->Fill();
  cout << "nbytes written:  " << nbytes << endl;
  reset();
}

void reset() {
  reset_pressures();
  reset_powers();
  //reset_temperatures();
}

void reset_pressures() {
  presCF4->Reset();
  presChamb->Reset();
}
void reset_temperatures() {
  //temp_tca_ptr->Clear();
}

void reset_powers() {
  powerScroll->Reset();
}

// open a log file and set the time for the next logfile to be opened
void open_log() {
  cout << "open_log()" << endl;
  //struct tm btime;  // broken down time

  // close any previously opened log file
  // (files should only be opened in the case of a rollover
  if (logfile != 0) {
    if (logfile->IsOpen()) {
      cout << "logfile is open.  closing...." << endl;
      logfile->Close();
    }
  } else {
    cout << "logfile is NOT open" << endl;
  }

  // compute the time for the next logfile rollover
  set_logfiledue();

  // make the filename for the log file
  // format is logYYYYMMDD.root
  // Data to be saved in /data/4shooter/log/ on rubin
  TTimeStamp *tst = new TTimeStamp();  // GMT time (don't use TDatime)
  logfilename  = "log";
  logfilename += tst->GetDate();
  logfilename += "_testing";
  logfilename += ".root";
  TString logfilename_full("");
  logfilename_full += logfiledir + logfilename;
  cout << "logfilename = " << logfilename << endl;
  cout << "logfilename_full = " << logfilename_full << endl;

  // open a ROOT file for logging.
  // If file already exists, append.  If not, then create.
  logfile = new TFile(logfilename_full, "UPDATE", "example log file");

  // check that file opened successfully
  if (verbose > 1) cout << "logfile->IsZombie() = " << logfile->IsZombie() << endl;
  if (logfile->IsZombie()) {    // open failed...
    cout << "logfile is zombie..." << endl;
    cout << "need to do something intelligent here..." << endl;
    exit(-1);
  }
  // a single tree holds all logged parameters
  logtree = new TTree("environment", "environment");

  // each parameter is logged to a separate branch
  logtree->Branch("presCF4",   "DmtpcLogEntry", &presCF4,   32000, SPLIT);
  logtree->Branch("presChamb", "DmtpcLogEntry", &presChamb, 32000, SPLIT);
  //presCF4_tca_ptr->BypassStreamer();  // took this from the $ROOTSYS/tutorials/tree/tcl.C example....
  //logtree->Branch("temperature", "TClonesArray", &temp_tca_ptr, 32000, SPLIT);
  logtree->Branch("powerScroll", "DmtpcLogEntry", &powerScroll, 32000, SPLIT);

}

void do_system_logging() {
  // handle opening logfile on startup and roll-over of log files
  if (verbose > 1) {
    cout << "do_system_logging()" << endl;
    cout << "dmtpctime, logfile_due = " << dmtpctime << ", " << logfile_due << endl;
  }
  if (dmtpctime >= logfile_due) {
    open_log();
  }

  //if (dmtpctime >= temp_due) {
  //  get_temperatures();
  //  temp_due += temp_interval;
  //}

  if (dmtpctime >= pres_due) {
    get_pressures();
    pres_due += pres_interval;
  }

  if (dmtpctime >= powers_due) {
    get_powers();
    powers_due += power_interval;
  }

  if (dmtpctime >= treefill_due) {
    fill_tree();
    treefill_due += treefill_interval;
  }

}

void do_waiting() {
  if (verbose > 1) cout << "do_waiting()" << endl;
  usleep(USLEEP_US);  // argument is microseconds
  //usleep(1000000);  // argument is microseconds
  set_dmtpctime();
}

/*---------------------------------------------------------------------------
init_time():  Initialize the time base and the current time

global out:  dmtpcbase, dmtpcctime are set to the current time
---------------------------------------------------------------------------*/
void init_time() {
  // initialize time
  struct timeval curtime;
  gettimeofday(&curtime, NULL);
  dmtpctbase = dmtpcctime = (ULong64_t)curtime.tv_sec*DMTPC_TPS
    + curtime.tv_usec/(1000000 / DMTPC_TPS);
}


/*---------------------------------------------------------------------------
init_loggable_params():  Initialize the variables to be logged
---------------------------------------------------------------------------*/
void init_loggable_params() {
  init_pressure_params();
}

void init_pressure_params() {
  
}

/*---------------------------------------------------------------------------
DetectorLogger(int nloop=0):  "infinite" loop for logging of system parameters

  if nloop = 0 then run indefinitely 
  set nloop to a positive value to execute that many iterations
  Each iteration takes at least USLEEP_US microseconds.

  By default, the logging code runs in an infinite loop.  To stop
  the logging code, simply create a file named .stoplog.  For
  example:

     > touch .stoplog 

  will stop the logger (as long as the file permissions are ok and you've done
  it in the correct directory).
---------------------------------------------------------------------------*/
int DetectorLogger(int nloop=0) {

  init_loggable_params();

  init_time();
  cout << "dmtpctbase = " << dmtpctbase << endl;

  int ii=0;
  reset();  // start fresh
  bool ok_to_log = true;
  while (ok_to_log) {
    if (verbose > 0) cout << "ii, dmtpctime = " << ii << ", " << dmtpctime << endl;
    do_system_logging();
    do_waiting();

    ii++;
    
    if ( (nloop > 0) and (ii >= nloop) ) { ok_to_log = false; }
    if (stop_requested()) { ok_to_log = false; }
    if (verbose > 1) cout << "       ok_to_log = " << ok_to_log << endl;
  }
  fill_tree();

  logtree->Print();
  logtree->Write();
  logfile->Close();

  cout << "logfilename = " << logfilename << endl;

  return EXIT_SUCCESS;
}

bool stop_requested() {
  // check to see if a file called .logstop exists
  TString stopfile(".stoplog");
  ifstream ifile(stopfile);
  if (ifile) return true;
  return false;
}

void read(TString filename) {
  cout << "open file" << endl;
  TFile *f = new TFile(filename);
  cout << "open tree" << endl;
  TTree *t = (TTree*)f->Get("environment");

  cout << "make tca" << endl;
  TClonesArray *tca = new TClonesArray("DmtpcLoggableParam", MAX_ENTRIES);
  TClonesArray *tca_presChamb = new TClonesArray("DmtpcLoggableParam", MAX_ENTRIES);
  cout << "branch address" << endl;
  t->GetBranch("presCF4")->SetAutoDelete(kFALSE);
  t->GetBranch("presChamb")->SetAutoDelete(kFALSE);
  t->SetBranchAddress("presCF4", &tca);
  t->SetBranchAddress("presChamb", &tca_presChamb);

  // loop over tree
  for (int ii=0; ii<t->GetEntries(); ii++) {
    tca->Clear();
    cout << "Event:  " << ii << endl;
    t->GetEvent(ii);
    // loop over tca for CF4 pressure
    for (int jj=0; jj<tca->GetEntries(); jj++) {
      cout << "  [" << tca->At(jj) << "]" << endl;
      cout << "  ";
      ((DmtpcLoggableParam*)tca->At(jj))->print();
    }
    // chamber pressure
    for (int jj=0; jj<tca_presChamb->GetEntries(); jj++) {
      cout << "  [" << tca_presChamb->At(jj) << "]" << endl;
      cout << "  ";
      ((DmtpcLoggableParam*)tca_presChamb->At(jj))->print();
    }
  }

  f->Close();
  delete tca;
}

// sample code to view logged data
TTree *readTree(TString filename) {

  TFile *f = new TFile(filename);
  TTree *t = (TTree*)f->Get("environment");

  //t->Draw("presCF4->GetValue():presCF4->GetTimeStamp()->GetTime()");
  return t;
}
 
void plotPres(TString filename) {

  TTree *tt = readTree(filename);
  //TCanvas *c = new TCanvas("junk", "junk");
  //c->Divide(2,1);
  //c->cd(1);
  tt->Draw("presCF4->GetValue():presCF4->GetTimeStamp()->GetTime()");
  //c->cd(2);
  //tt->Draw("presChamb->GetValue():presChamb->GetTimeStamp()->GetTime()");
  
}

void readScroll(TString filename) {
  TTree *tt = readTree(filename);
  cout << "make tca" << endl;
  TClonesArray *tca = new TClonesArray("DmtpcLoggableParam", MAX_ENTRIES);
  cout << "branch address" << endl;
  tt->GetBranch("powerScroll")->SetAutoDelete(kFALSE);
  tt->SetBranchAddress("powerScroll", &tca);

  // loop over tree
  for (int ii=0; ii<tt->GetEntries(); ii++) {
    tca->Clear();
    cout << "Event:  " << ii << endl;
    tt->GetEvent(ii);
    // loop over tca 
    cout << "tca->GetEntries() = " << tca->GetEntries() << endl;
    for (int jj=0; jj<tca->GetEntries(); jj++) {
      cout << "  [" << tca->At(jj) << "]" << endl;
      cout << "  ";
      ((DmtpcLoggableParam*)tca->At(jj))->print();
    }
  }

}

void graphPresCF4(TString filename) {
  TFile f(filename);
  TTree *t = (TTree*)f.Get("environment");

  TGraph *gPresCF4 = new TGraph();

  DmtpcLogEntry *dle = new DmtpcLogEntry();
  cout << "setbranchaddress" << endl;
  t->SetBranchAddress("presCF4", &dle);
  Int_t nPts = 0;
  for (int ii=0; ii<t->GetEntries(); ii++) {
    t->GetEvent(ii);
    cout << dle->GetNSamples() << endl;
    dle->print();
    for (Int_t jj=0; jj<dle->GetNSamples(); jj++) {
      gPresCF4->SetPoint(nPts, 
			 dle->GetTimeStamp(jj)->GetTime(), 
			 dle->GetValue(jj));
      nPts++;
    }
  }
  gPresCF4->Draw("alp");
}
void graphPresChamb(TString filename) {
  TFile f(filename);
  TTree *t = (TTree*)f.Get("environment");

  TGraph *gPresChamb = new TGraph();

  DmtpcLogEntry *dle = new DmtpcLogEntry();
  cout << "setbranchaddress" << endl;
  t->SetBranchAddress("presChamb", &dle);
  Int_t nPts = 0;
  for (int ii=0; ii<t->GetEntries(); ii++) {
    t->GetEvent(ii);
    cout << dle->GetNSamples() << endl;
    dle->print();
    for (Int_t jj=0; jj<dle->GetNSamples(); jj++) {
      gPresChamb->SetPoint(nPts, 
			 dle->GetTimeStamp(jj)->GetTime(), 
			 dle->GetValue(jj));
      nPts++;
    }
  }
  gPresChamb->Draw("alp");
}


void read_powers() {
  // run the python script
  system("python getpower.py");
  // open the output file: power.tmp
  ifstream fin("power.tmp");
  TString portname, portstatus;

  //  1. scroll pump
  //  2. ccd0
  //  3. ccd1
  //  4. ccd2
  //  5. ccd3
  //  6. AC/DC 0 --> 4.5" Gate Valve
  //  7. AC/DC 1 --> 2.75" Gate Valve
  //  8. AC/DC 2 --> Gas input valve
  //  9. AC/DC 3 --> Scroll valve
  // 10. 
  // 11. 
  // 12. 
  // 13. 
  // 14. 
  // 15. 
  // 16. Mass Flow Controller
  
  while (!fin.eof()) {
    fin >> portname >> portstatus;
    cout << "portname, status = " << portname << ", " << portstatus << endl;
    if (portname == "Port1") { }
    else if (portname == "Port2") {}
    else if (portname == "Port3") {}
    else if (portname == "Port4") {}
    else if (portname == "Port5") {}
    else if (portname == "Port6") {}
    else if (portname == "Port7") {}
    else if (portname == "Port8") {}
    else if (portname == "Port9") {}
    else if (portname == "Port10") {}
    else if (portname == "Port11") {}
    else if (portname == "Port12") {powerScroll_param->SetValue(valOfPowerStatus(portstatus)); }
    else if (portname == "Port13") {}
    else if (portname == "Port14") {}
    else if (portname == "Port15") {}
    else if (portname == "Port16") {}
  }
  cout << "powerScroll_param->GetValue() = " << powerScroll_param->GetValue() << endl;
}

Float_t valOfPowerStatus(TString status) {
  if (status == "Off") { return 0.0; }
  else if (status == "On") { return 1.0; }
  else { return -1.0; }
}
