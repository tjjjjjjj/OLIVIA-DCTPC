#include <iostream>
using std::cout;
using std::endl;

#include "TROOT.h"
#include "TSystem.h"
#include "TSocket.h"
#include "MaxCamChannel.hh"
#include <cmath>

char * default_host = "18.146.0.100"; 

// to do list:
// 1. flag unset pressure, and check that something reasonable happens
// 2  make sure pressure changes when pump is on and you open the valve
// 3  " for gas
// 4  find the right value of dP
// 5  correct the exit logic
// 6  figure out what the running script does: "more MaxCam/script/setpress"

// Pressure tolerance (Torr)
double dP = 0.4; 
// Pump down pressure (Torr)
double DP = 0.02;

int main(int argc, char* argv[]) {
  
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

  // notify everyone that you're doing something
  redsock->SendRaw("pset 8 1\n", 10);
  gSystem->Sleep(500);
  cout<<"debug: check if socket 8 is open now ... "<<endl;

  // close all valves at the begining
  bool isGasOpen=false;
  if (debug) {
    redsock->SendRaw("pset 7 1\n", 10); // test
    cout<<"debug: check if socket 7 is open now ... "<<endl;
    getchar();
  } else {
    redsock->SendRaw("pset 4 0\n", 10); // close gas valve
    gSystem->Sleep(500);
  }

  bool isPumpOpen=false;
  if (debug) {
    redsock->SendRaw("pset 7 0\n", 10); // test
    cout<<"debug: check if socket 7 is closed now ... "<<endl;
    getchar();
  } else {
    redsock->SendRaw("pset 6 0\n", 10); // close pump valve
    gSystem->Sleep(500);
  }

  // We always want to open the valve between chamber and bottle valve at the beginning
  bool isGasChamberOpen=true;
  if (debug) {
    redsock->SendRaw("pset 7 1\n", 10); // test
    cout<<"debug: check if socket 7 is open now ... "<<endl;
    getchar();
  } else {
    redsock->SendRaw("pset 5 1\n", 10); // close gas valve
    gSystem->Sleep(500);
  }

  redsock->SendRaw("\r", 2); // clear
  gSystem->Sleep(500);

  cout<<"Start loop"<<endl;
  int loop_counter2 = 0;
  MaxCamChannel ch("pressure");
  ch.setBusy(true); 
  while(1){
    //read presure from database
    ch.readFromDB("pressure");
    double currentP = ch.currentValue;
    double setP = ch.setValue;

    if (adjust) {
      break;
    }

    if (currentP<DP){
      break;
    }

    if (currentP>DP && !isPumpOpen){
      //cout<<"Opening pump valve"<<endl;
      redsock->SendRaw("pset 4 0\n", 10); //make sure gas valve is closed
      gSystem->Sleep(500);
      redsock->SendRaw("pset 6 1\n", 10); //open pump valve
      gSystem->Sleep(500);
      isGasOpen=false;
      isPumpOpen=true;
    }
  
    else if (currentP>DP && isPumpOpen){
      cout<<"    Pumping . . ."<<endl;
      redsock->SendRaw("pset 4 0\n",10); //make sure gas valve is closed
      gSystem->Sleep(500);
      isGasOpen=false;
      isPumpOpen=true;
    }
    gSystem->Sleep(3500);
    loop_counter2++;
    
  }

  cout << "Start loop"<<endl;
  int loop_counter = 0;
  while (1) {

    // read pressure from data base
    ch.readFromDB("pressure");    
    double currentP = ch.currentValue;
    double setP = ch.setValue;

    cout << "Pressure ("<<loop_counter<<") : current="<<currentP 
	 << "  set=" << setP << "   dP=" << dP 
	 << "   valves: gas="<<isGasOpen<<" pump="<<isPumpOpen
	 << endl;
    

    //If pressure is OK, make sure all valves are closed 
    if (fabs(currentP-setP)<dP) {
      cout<<"Pressure OK"<<endl;
      redsock->SendRaw("pset 6 0\n", 10);
      gSystem->Sleep(500);
      redsock->SendRaw("pset 4 0\n", 10);
      gSystem->Sleep(500);
      redsock->SendRaw("pset 5 0\n", 10);
      cout<<"re-checking pressure (2 minutes)"<<endl;
      gSystem->Sleep(120000);
      ch.readFromDB("pressure");
      double currentP = ch.currentValue;
      double setP = ch.setValue; 
      if (fabs(currentP-setP)<dP) {
	cout<<"Pressure OK"<<endl;
      break;
      } else {
	cout<<"re-adjusting pressure"<<endl;
	//reopen chamber gas
	redsock->SendRaw("pset 5 1\n",10);
	gSystem->Sleep(500);
      }
    }

    else if ( setP-currentP>dP) {
      cout << "Open gas valve" << endl;
      if (debug) {
	redsock->SendRaw("pset 7 0\n", 10); // test
	cout<<"check if socket 7 is closed now ... "<<endl;
	getchar();
	gSystem->Sleep(500);
      } else {
	cout<<"     setP-currentP = "<<setP-currentP<<endl;
	redsock->SendRaw("pset 6 0\n", 10); // make sure pump is closed
	gSystem->Sleep(500);
	redsock->SendRaw("pset 4 1\n", 10); // open gas valve 
	if ( setP-currentP>10) {gSystem->Sleep(8000);
	} else {
	gSystem->Sleep(500);
	}
	redsock->SendRaw("pset 4 0\n", 10); //close gas valve
	cout<<"     allowing database to update . . ."<<endl;
	gSystem->Sleep(2500);
      }
      isPumpOpen=false;
      isGasOpen=false;
    }

    else if ( currentP-setP>dP) {
      cout << "Open pump valve" << endl;
      if (debug) {
	redsock->SendRaw("pset 7 0\n", 10); // test
	cout<<"check if socket 7 is closed now (3) ..."<<endl;
	getchar();
	gSystem->Sleep(500);
      } else {
	cout<<"     currentP-setP = "<<currentP-setP<<endl;
	redsock->SendRaw("pset 4 0\n", 10); // make sure gas is closed
	gSystem->Sleep(500);
	redsock->SendRaw("pset 6 1\n", 10); // open pump valve
	gSystem->Sleep(1500);
	redsock->SendRaw("pset 6 0\n",10); // close pump valve
	cout<<"     allowing database to update . . ."<<endl;  
	gSystem->Sleep(2500);
      }
      isGasOpen=false;
      isPumpOpen=false;
    }
    gSystem->Sleep(2500);
    loop_counter++;
  }

  // finish up, close all values
  redsock->SendRaw("pset 4 0\n", 10); // make sure gas is closed
  gSystem->Sleep(500);
  redsock->SendRaw("pset 6 0\n", 10); // make sure pump valve is closed
  gSystem->Sleep(500);
  redsock->SendRaw("pset 5 0\n",10); //close chamber gas valve

  // notify everyone that you're done doing something
  redsock->SendRaw("pset 8 0\n", 10); 
  gSystem->Sleep(500);
  cout<<"check if sockets 4, 6, 8 are closed now ... "<<endl;

  cout << "Closing socket"<<endl;
  redsock->Close();
  delete redsock;
  ch.setBusy(false); 
  return 0;
}
