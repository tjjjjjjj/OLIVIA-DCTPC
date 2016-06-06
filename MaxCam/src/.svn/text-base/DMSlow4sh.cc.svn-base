#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::setprecision;

#include "TROOT.h"
#include "TSystem.h" 
#include "MaxCamChannel.hh"
#include <cmath>

void readStandardChannel(MaxCamChannel *ch, 
			 Double_t level_of_change_that_is_interesting,
			 Double_t scale_factor=1.0,
			 Bool_t ForceRead=true);

MaxCamChannel * pressure_cdg; 
MaxCamChannel * pressure_bpg; 
MaxCamChannel * pressure_convectron; 

MaxCamChannel * hvstatus;
MaxCamChannel * wire_i;
MaxCamChannel * wire_hv;
MaxCamChannel * mesh_i;
MaxCamChannel * mesh_hv;

MaxCamChannel * temp0;
MaxCamChannel * temp1;

// the mass flow controller
MaxCamChannel * mfc;

// the caen readout and control channels
MaxCamChannel * caen_n1471_hv0;
MaxCamChannel * caen_n1471_hv1;
MaxCamChannel * caen_n1471_hv2;
MaxCamChannel * caen_n1471_hv3;

Int_t debug;

void readMFC(){
  if(debug) cout << "in readMFC() ..." << endl;

  double old_mfc_val = mfc->currentValue;
  // have to recall the setval from the databases' last entry
  mfc->readValueFromDB("mfc",2); 
  // read the current value
  mfc->readADC(100);

  if(debug) cout << mfc->GetName() 
		 << ": " 
		 << mfc->currentValue   
		 << "  +/- " 
		 << mfc->currentValueRMS   
		 <<  " set=" 
		 << mfc->setValue 
		 << " fabs(old_mfc_val - mfc->currentValue)=" 
		 << fabs(old_mfc_val - mfc->currentValue)
		 << " old_mfc_val="
		 << old_mfc_val
		 << endl; 
  
  if (fabs(old_mfc_val - mfc->currentValue) > 0.005) mfc->saveValueToDB();

  // set the mfc to the current set-point
  mfc->setDAC();

  if(debug) cout << "leaving readMFC() ..." << endl;
}

void readConvectron() {

  // convectron gauge (RANGE?)
  double old_convectron_val = pressure_convectron->currentValue;
  // have to recal the setval from the databases' last entry
  pressure_convectron->readValueFromDB("pressure_convectron",2); 
  // read the current value
  pressure_convectron->readConvectron(NULL,false);

  // for debugging the pressure
  if(debug) cout << "old_convectron_val=" 
		 << setprecision(6) << old_convectron_val 
		 << " pressure_convectron->currentValue=" 
		 << setprecision(6) << pressure_convectron->currentValue 
		 << " delta_pressure_convectron=" 
		 << setprecision(6) << fabs(old_convectron_val-pressure_convectron->currentValue) 
		 << endl;

  // looks up the name if one isn't given, so this is ok
  if (fabs(old_convectron_val - pressure_convectron->currentValue) > 0.001) pressure_convectron->saveValueToDB();
}

void readCDG() {

  // capacitance diaphragm gauge (0-200 torr, accurate)
  double old_cdg_val = pressure_cdg->currentValue;
  // have to recal the setval from the databases' last entry
  pressure_cdg->readValueFromDB("pressure_cdg",2); 
  // read the current value
  pressure_cdg->readInficonController(NULL,2,9600);
  
  // for debugging the pressure
  if(debug) cout << "old_cdg_val=" 
		 << setprecision(6) << old_cdg_val 
		 << " pressure_cdg->currentValue=" 
		 << setprecision(6) << pressure_cdg->currentValue 
		 << " delta_pressure_cdg=" 
		 << setprecision(6) << fabs(old_cdg_val-pressure_cdg->currentValue) 
		 << endl;
  
  // looks up the name if one isn't given, so this is ok
  if (fabs(old_cdg_val - pressure_cdg->currentValue) > 0.002) pressure_cdg->saveValueToDB();

}

void readBPG() {

  // (dual gauge, atmosphere to 1e-11 Torr)
  double old_bpg_val = pressure_bpg->currentValue;
  // have to recal the setval from the databases' last entry
  pressure_bpg->readValueFromDB("pressure_bpg",2); 
  // read the current value
  pressure_bpg->readInficonController(NULL,1,9600);

  // for debugging the pressure
  if(debug) cout << "old_bpg_val=" 
		 << setprecision(6) << old_bpg_val 
		 << " pressure_bpg->currentValue=" 
		 << setprecision(6) << pressure_bpg->currentValue 
		 << " delta_pressure_bpg=" 
		 << setprecision(6) << fabs(old_bpg_val-pressure_bpg->currentValue) 
		 << endl;

  if( old_bpg_val!=0 ){
    cout << "old_bpg_val=" << old_bpg_val << endl;
    cout << "pressure_bpg->currentValue=" << pressure_bpg->currentValue << endl;
    cout << "fabs(old_bpg_val - pressure_bpg->currentValue)/fabs(old_bpg_val)="
	 << fabs(old_bpg_val - pressure_bpg->currentValue)/fabs(old_bpg_val)
	 << endl;
    if (fabs(old_bpg_val - pressure_bpg->currentValue)/fabs(old_bpg_val) > 0.005) pressure_bpg->saveValueToDB();
  } else {
    if (fabs(old_bpg_val - pressure_bpg->currentValue) > 0.005) pressure_bpg->saveValueToDB();
  }
}

void readWithoutRamping() {
  double old_val = wire_hv->currentValue;
  wire_hv->readADC();
  if(debug) cout << wire_hv->GetName() << ": " << wire_hv->currentValue   << "  +/- " << wire_hv->currentValueRMS   <<  " set=" << wire_hv->setValue << endl; 
  
  wire_hv->setValue = 0;
  wire_hv->currentValueRMS = 0;
  
  if (fabs(old_val - wire_hv->currentValue) > 0.01) wire_hv->saveValueToDB();

  old_val = mesh_hv->currentValue;
  mesh_hv->readADC();
  if(debug) cout << mesh_hv->GetName() << ": " << mesh_hv->currentValue   << "  +/- " << mesh_hv->currentValueRMS   <<  " set=" << mesh_hv->setValue << endl; 
  
  mesh_hv->setValue = 0;
  mesh_hv->currentValueRMS = 0;
  
  if (fabs(old_val - mesh_hv->currentValue) > 0.01) mesh_hv->saveValueToDB();
}

void rampFromDB() {

  double deltaHV=0.005; // 5V
  
  // change vol if set values different
  // or, update current value in database

  hvstatus->readFromDB("hvstatus");
  bool hvON = hvstatus->setValue>0 ? true : false;

  if(debug) cout << "HV STATUS: " << hvON << endl; 

  if (!hvON) {
    wire_hv->ramp(0.0);
    double last_wire_hv_val = wire_hv->currentValue;
    wire_hv->readValueFromDB("wire_hv",2); 
    wire_hv->currentValue = last_wire_hv_val; 
    wire_hv->saveValueToDB("wire_hv"); 
    wire_hv->print();

    mesh_hv->ramp(0.0);
    double last_mesh_hv_val = mesh_hv->currentValue;
    mesh_hv->readValueFromDB("mesh_hv",2); 
    mesh_hv->currentValue = last_mesh_hv_val; 
    mesh_hv->saveValueToDB("mesh_hv"); 
    mesh_hv->print();

    

    return;
  } 

  double wire_hv_old_setval = wire_hv->readValueFromDB("wire_hv",2);
  double wire_hv_old_val = wire_hv->readValueFromDB("wire_hv",0);
  
  wire_hv->readADC(1000);
  if (fabs(wire_hv_old_setval-wire_hv->currentValue)>1e-4)      { 
    wire_hv->rampFromDB("wire_hv"); 
  }
  else if (fabs(wire_hv_old_val-wire_hv->currentValue)>deltaHV) { 
    wire_hv->saveValueToDB("wire_hv"); 
  } 
  gSystem->Sleep(500);
  
  if(debug) cout << "wire_hv->setValue=" << wire_hv->setValue << " wire_hv->currentValue=" << wire_hv->currentValue << endl;
  
  double mesh_hv_old_setval = mesh_hv->readValueFromDB("mesh_hv",2);
  double mesh_hv_old_val = mesh_hv->readValueFromDB("mesh_hv",0);
  mesh_hv->readADC(1000);
  if (fabs(mesh_hv_old_setval-mesh_hv->currentValue)>1e-4)      { 
    mesh_hv->rampFromDB("mesh_hv"); 
  } 
  else if (fabs(mesh_hv_old_val-mesh_hv->currentValue)>deltaHV) { 
    mesh_hv->saveValueToDB("mesh_hv"); 
  } 

  gSystem->Sleep(500);
}

void readStandardChannel(MaxCamChannel *ch, 
			 Double_t level_of_change_that_is_interesting,
			 Double_t scale_factor,
			 Bool_t ForceRead){

  double old_val = ch->currentValue;
  // have to recal the setval from the databases' last entry
  ch->readValueFromDB(ch->GetName(),2); 

  ch->readADC();
  ch->currentValue*=scale_factor;
  ch->currentValueRMS*=scale_factor;

  cout << ch->GetName() << ": " << ch->currentValue   << "  +/- " << ch->currentValueRMS   <<  " set=" << ch->setValue << endl; 

  // implement this, but later
  ch->currentValueRMS = 0;
  
  if ( ( fabs(old_val - ch->currentValue) > level_of_change_that_is_interesting ) ||
       ForceRead ) ch->saveValueToDB();

}

int main(int argc, char *argv[]) {

  debug=false;

  // initialize the pressure channels
  pressure_cdg = new MaxCamChannel("pressure_cdg"); 
  pressure_bpg = new MaxCamChannel("pressure_bpg"); 
  pressure_convectron = new MaxCamChannel("pressure_convectron"); 

  // initialize the high-voltage related channels
  wire_i = new MaxCamChannel("wire_i");
  wire_hv = new MaxCamChannel("wire_hv");
  mesh_i = new MaxCamChannel("mesh_i");
  mesh_hv = new MaxCamChannel("mesh_hv");
  hvstatus = new MaxCamChannel("hvstatus","",-1,-1,"/etc/dbaccess.txt");

  // initialize the caen-related channels
  caen_n1471_hv0 = new MaxCamChannel("caen_n1471_hv0");
  caen_n1471_hv1 = new MaxCamChannel("caen_n1471_hv1");
  caen_n1471_hv2 = new MaxCamChannel("caen_n1471_hv2");
  caen_n1471_hv3 = new MaxCamChannel("caen_n1471_hv3");

  // initialize the temperature channels
  temp0 = new MaxCamChannel("temp0");
  temp1 = new MaxCamChannel("temp1");

  //    mass flow controller
  mfc = new MaxCamChannel("mfc");

  // default sleep time
  int nsleep=500;
  // really short sleep time for debugging purposes
  //  int nsleep=1;
  int nskip=1; // read voltage every nsleep*nskip msec
  if (argc>1) {
    nsleep=atoi( argv[1] );
    if(debug) cout << argc << "   " << (argv[1]) << endl;
  }
  if(debug)  cout << "Sleep set to " << nsleep << " ms" << endl;
  
  int nread=0;
  while (1) {
    // read out the pressure gauges, if they're not busy
    if( !( pressure_convectron->isBusy() ) ); readConvectron();
    if( !( pressure_cdg->isBusy() ) ) readCDG();
    if( !( pressure_bpg->isBusy() ) ) readBPG();

    if( !( mfc->isBusy() ) ) readMFC();

    if (nread++ % nskip == 0) {
      // read out the HV currents
      readStandardChannel(wire_i,0.01);
      readStandardChannel(mesh_i,0.01);

      // should be 0.01
      readStandardChannel(caen_n1471_hv0,0.1);
      gSystem->Sleep(10);
      readStandardChannel(caen_n1471_hv1,0.1);
      gSystem->Sleep(10);
      readStandardChannel(caen_n1471_hv2,0.1);
      gSystem->Sleep(10);
      readStandardChannel(caen_n1471_hv3,0.1);
      gSystem->Sleep(10);

      caen_n1471_hv0->setDAC();
      caen_n1471_hv1->setDAC();
      caen_n1471_hv2->setDAC();
      caen_n1471_hv3->setDAC();

      //      readWithoutRamping();
      if(debug) cout << "right before rampFromDB()..." << endl;
      rampFromDB();
      if(debug) cout << "right after rampFromDB()..." << endl;

      readStandardChannel(temp0,0.01,100.);
      readStandardChannel(temp1,0.01,100.);
    }

    gSystem->Sleep(nsleep);
  }

  return 0;
}


