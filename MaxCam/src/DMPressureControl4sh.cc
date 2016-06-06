#include <iostream>
using std::cout;
using std::endl;

#include "TROOT.h"
#include "TSystem.h"
#include "TSocket.h"
#include "MaxCamChannel.hh"
#include <cmath>

// the ip of the synaccess; cannon.lns.mit.edu for
// the 4-shooter
char * default_host = "cannon.lns.mit.edu"; 

// Pressure tolerance (Torr)
double dP = 0.9; 
// Pump down pressure (Torr)
double pump_down_minimum_pressure = 0.01;

// pressure channels
MaxCamChannel * pressure_cdg;
MaxCamChannel * pressure_bpg;
MaxCamChannel * pressure_convectron;

// the mass flow controller
MaxCamChannel * mfc;

// to do
// 1) incorporate busy status stuff into the fill controls

void readPressure();
void savePressure();
void readMFC();

void readPressure() {
  cout << "in readPressure() ..." << endl;

  cout << "here1" << endl;
  // this is for reading the cdg gauge
  // have to recall the setval from the databases' last entry
  pressure_cdg->readValueFromDB("pressure_cdg",2); 
  // read the current value
  cout << "here2" << endl;
  pressure_cdg->readInficonController(NULL,2,9600);
  cout << "here3" << endl;
  
  cout << "pressure_cdg->currentValue=" << pressure_cdg->currentValue << " pressure_cdg->setValue=" << pressure_cdg->setValue << endl;
  cout << pressure_cdg->GetName() << ": " << pressure_cdg->currentValue   << "  +/- " << pressure_cdg->currentValueRMS   <<  " set=" << pressure_cdg->setValue << endl;
  
  cout << "here4" << endl;

  cout << "leaving readPressure() ..." << endl;
}

void savePressure() {
  cout << "in savePressure() ..." << endl;
  pressure_cdg->saveValueToDB();
}


void readMFC(){
  cout << "in readMFC() ..." << endl;

  // have to recall the setval from the databases' last entry
  mfc->readValueFromDB("mfc",2); 
  // read the current value
  mfc->readADC(100);

  cout << "mfc->currentValue=" << mfc->currentValue << " mfc->setValue=" << mfc->setValue << endl;
  cout << mfc->GetName() << ": " << mfc->currentValue   << "  +/- " << mfc->currentValueRMS   <<  " set=" << mfc->setValue << endl; 

  cout << "leaving readMFC() ..." << endl;
}

void saveMFC(){
  cout << "in saveMFC() ..." << endl;
  // probably don't actually want this here...
  //  readMFC();
  mfc->saveValueToDB();
}

void setMFC(){

  // read the current value right before we change the set value
  mfc->readADC(100);
  // set the flow rate of the mfc to whatever's currently in the setval
  mfc->setDAC();
  // save the set value to the database
  saveMFC();
}

int main(int argc, char* argv[]) {

  cout << "In DMPressureControl4sh.cc main()..." << endl;

  // initialize channels 
  //    pressure channels
  pressure_cdg = new MaxCamChannel("pressure_cdg");
  pressure_bpg = new MaxCamChannel("pressure_bpg");
  pressure_convectron = new MaxCamChannel("pressure_convectron");

  //    mass flow controller
  mfc = new MaxCamChannel("mfc");

  // set the pressure channels to busy, so that we don't try 
  // to dual read the pressure while filling; this leads to 
  // freeze-outs otherwise when two devices try to access 
  // the same pressure gauge at the same time over serial
  pressure_cdg->setBusy(true); 
  pressure_bpg->setBusy(true);
  pressure_convectron->setBusy(true);
  mfc->setBusy(true);

  // need to give the slow control a little time to realize the busy 
  // flags have been set
  gSystem->Sleep(500);
  
  bool debug = false;
  bool adjust = false;

  if (argc>1) {
    TString adjust_flag(argv[1]);
    if (adjust_flag=="adjust") {
      adjust = true;
      cout<<"Adjusting pressure . . ." <<endl;
    }
  } else {
    adjust = false;
  }  

  if (argc>1) {
    TString debug_flag(argv[1]);
    cout<<"program arguments: 1 "<<argv[1]<<endl;
    cout<<"debug flag: "<<debug_flag<<endl;
    if (debug_flag=="debug") {
      debug = true;
      cout<<"***** running in debug mode! *****"<<endl;
    } 
  } else {
    debug = false;
  }

  cout << "Open socket" << endl;
  //TSocket *redsock = new TSocket("powmgmt.lns.mit.edu", 23);
  char * host; 
  if (argc>2) host = argv[2];
  else host = default_host; 
  TSocket *redsock = new TSocket(host, 23);
  if (!redsock->IsValid()) {
      cout << "Socket not valid " << endl;
      return 1;
  } else {
    if (debug) {cout<<"Successfully opened socket 18.146.0.100"<<endl;}
  }

  // check initial message
  char rmsg[2048];
  cout << "Check initial message"<<endl;
  int nbytes = redsock->RecvRaw( rmsg, 2048, kDontBlock);
  cout << rmsg << endl;
  nbytes = redsock->RecvRaw( rmsg, 2048, kDontBlock);
  cout << rmsg << endl;

  // keep track of whether or not, in the refill section of the code,
  // we've told the rough pump to be on, or off.
  Bool_t isRoughPumpOn=false;
  cout<<"Start loop"<<endl;
  int loop_counter = 0;
  //  ch.setBusy(true); 
  /*
  while(1){

    if (adjust) {
      // we're just trying to match the set value, don't want 
      // to do any pumping
      break;
    }
    
    //read the current pressure
    readCV1Pressure();
    double currentP = pressure_cv1->currentValue;
    double setP = pressure_cv1->setValue;

    cout << "Pressure ("<<loop_counter<<") : current="<<currentP 
	 << "  set=" << setP << "   dP=" << dP 
	 << endl;
    
    cout << "currentP=" << currentP << " setP=" << setP << endl;
    
    cout << "pump_down_minimum_pressure=" << pump_down_minimum_pressure << endl;

    if(currentP<pump_down_minimum_pressure){

      cout << "(currentP<pump_down_minimum_pressure)" << endl;
      
      redsock->SendRaw("pset 14 0\n", strlen("pset 14 0\n")+1);
      // rough pump should be off, so 
      isRoughPumpOn=false;
      // we've reached our target, so break!
      break;

    } else {

      cout << "!(currentP<pump_down_minimum_pressure)" << endl;
      cout << "isRoughPumpOn=" << isRoughPumpOn << endl;

      // if the pump isn't on; turn it on.
      if(!isRoughPumpOn){
	cout << "... in if(!isRoughPumpOn) { ... }" << endl;
	redsock->SendRaw("pset 14 1\n", strlen("pset 14 1\n")+1);
      }

    }

    // amount to wait in between checking to see if we've
    // pumped down lower than the minimum pressure
    // 3500...
    gSystem->Sleep(500);
    ++loop_counter;
  }
  */
  
  // we're done with anything involving the pumps, so make sure that
  // the pump is off.  Send it the close command, just to be super redundant.
  //  redsock->SendRaw("pset 14 0\n", strlen("pset 14 0\n")+1);
  
  // The ``adjust'' pressure loop
  //  re-initialize loop_counter for use in the adjust loop
  loop_counter=0;
  cout << "Start the adjust loop"<<endl;
  while (1) {
    
    // read the current pressure
    readPressure();
    double currentP = pressure_cdg->currentValue;
    double setP = pressure_cdg->setValue;
    savePressure();
    
    cout << "Pressure ("<<loop_counter<<") : current="<<currentP 
	 << "  set=" << setP << "   dP=" << dP 
	 << endl;
    
    // if pressure is lower than the set pressure, open the mfc to full
    // first read the MFC to see what it's current flow rate is
    readMFC();
    double currentFlow = mfc->currentValue;
    double setFlow = mfc->setValue;
    // we just read the mfc, save the read value
    saveMFC();
    
    cout << "Flow ("<<loop_counter<<") : current="<<currentFlow
	 << "  set=" << setFlow
	 << endl;
    
    // so now we know the current flow rate and the current pressure.
    // fill until we edge just past the pressure set-point
    if(currentP>setP){
      cout << "(currentP>setP)" << endl;
      // stop filling by setting the flow rate to zero
      mfc->setValue=0.;
      setMFC();
      // we're done adjusting the pressure; break out of this loop
      break;
    } else {
      // we're not to the set-pressure yet.  Keep flowing.
      mfc->setValue=5.;
      setMFC();
    }
    
    gSystem->Sleep(100);
    ++loop_counter;
    // for debugging, can break here...
    //    if(loop_counter>20) break;
  }

  cout << "Closing socket"<<endl;
  redsock->Close();
  delete redsock;

  // we're done with all fill-related activities requested by the 
  // user; free up the pressure channels so they can be safely used
  // by other things, like the slow control
  pressure_cdg->setBusy(false); 
  pressure_bpg->setBusy(false);
  pressure_convectron->setBusy(false);
  mfc->setBusy(false);
  
  cout << "Leaving DMPressureControl4sh.cc main()..." << endl;
  return 0;
}
