#include <iostream>
#include <iomanip>
#include <cstdlib>
using std::cout;
using std::endl;
using std::setprecision;

#include "TROOT.h"
#include "TSystem.h" 
#include "TMath.h"
#include "MaxCamChannel.hh"
#include <cmath>

void readTemperature(MaxCamChannel *ch);
void readMKS974B(MaxCamChannel *ch, Double_t level_of_change_that_is_interesting, Bool_t ForceRead=true);

void readStandardChannel(MaxCamChannel *ch, 
			 Double_t level_of_change_that_is_interesting,
			 Double_t scale_factor=1.0,
			 Bool_t ForceRead=true);

MaxCamChannel * pressure_mks974B;
MaxCamChannel * temp0;

Int_t debug;

void readMKS974B(MaxCamChannel *ch, Double_t level_of_change_that_is_interesting, Bool_t ForceRead) {
  double old_val = ch->currentValue;
  ch->readValueFromDB(ch->GetName(),2);
  ch->readADC();

  // convert voltage to pressure
  ch->currentValue    = TMath::Power(10, 2.0*ch->currentValue-11.);
  ch->currentValueRMS = 0.0;

  cout << ch->GetName() << ": " << ch->currentValue   << "  +/- " << ch->currentValueRMS   <<  " set=" << ch->setValue << endl; 

  if ( ( fabs(old_val - ch->currentValue) > level_of_change_that_is_interesting ) ||
       ForceRead ) ch->saveValueToDB();
}

void readTemperature(MaxCamChannel *ch) {
  readStandardChannel(ch, 0.01, 100.);
}


void readStandardChannel(MaxCamChannel *ch, 
			 Double_t level_of_change_that_is_interesting,
			 Double_t scale_factor,
			 Bool_t ForceRead){

  double old_val = ch->currentValue;

  // have to recal the setval from the databases' last entry
  // this updates ch->setValue)
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

  debug=true;

  // initialize the pressure channels
  pressure_mks974B = new MaxCamChannel("pressure_mks974B");

  // initialize the temperature channels
  temp0 = new MaxCamChannel("temp0");

  // default sleep time (ms)
  int nsleep=5000;
  // really short sleep time for debugging purposes
  //  int nsleep=1;

  //int nskip=1; // read voltage every nsleep*nskip msec

  if (argc>1) {
    nsleep=atoi( argv[1] );
    if(debug) cout << argc << "   " << (argv[1]) << endl;
  }
  if(debug)  cout << "Sleep set to " << nsleep << " ms" << endl;
  
  int nread=0;
  while (1) {
    // read out the pressure gauges, if they're not busy
    //if( !( pressure_convectron->isBusy() ) ); readConvectron();

    readMKS974B(pressure_mks974B, 0.01);
    readTemperature(temp0);

    gSystem->Sleep(nsleep);
  }

  return 0;
}


