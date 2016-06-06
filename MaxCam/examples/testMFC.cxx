// check what's happening on a port viz-a-viz
// chmod a+rw /dev/ttyS*
// sudo chmod a+rw /dev/ttyS*
// ls -la /dev/ttyS*
// cat /dev/ttyS1
// cat /dev/ttyS0
// cat /dev/ttyS2

#include <unistd.h>
#include <termios.h>
#include <string>
#include <iostream>
#include <fcntl.h>

#include "../MaxCamMFC.hh"
#include "../MaxCamSerial.hh"
using namespace std;

void testMFC(void){
  // on the four shooter daq, there are multiple serial ports
  // so we got to make sure that we're pointing to the right 
  // one.  Can check which one is the right one through the
  // commands at the top of this code.
  MaxCamMFC *mfc = new MaxCamMFC("mfc","mfc","/dev/ttyS2");
  // we've opened up a connection to the MFC

  // tell controller to take all commands from serial
  //  mfc->setSource(1);
  // set units; 
  //  mfc->setUnits(23); //23==g/min (max flow rate is 0.494 g/min)
  // set the gas;
  //  mfc->setGas(7); // 7==90% Argon/10% Methane
  // zero the set-point 
  //  mfc->setFlowPoint(0.0);
  // open valve
  //  mfc->openValve();
  
  //  sleep(1);
  
  //  mfc->closeValve();

  // an example of a fill cycle
  /*
  sleep(1);

  // fill at 0.1 g/min for 2 minutes
  float flowrate=0.05; // g/min
  int flowtime=60*2; // seconds
  mfc->Fill(flowrate,flowtime);
  
  // let's leave the valve closed, to be safe
  mfc->closeValve();
  */
  
  //  mfc->setDebug(true);
  mfc->setPause(0.3);

  // all of the things you can read from the MFC
  cout << "mfc->getVersion()=" << mfc->getVersion() << endl;
  cout << "mfc->getFactor()=" << mfc->getFactor() << endl;
  cout << "mfc->getSerial()=" << mfc->getSerial() << endl;
  cout << "mfc->getFlowPoint()="<< mfc->getFlowPoint() << endl;
  cout << "mfc->getSource()=" << mfc->getSource() << endl;
  cout << "mfc->getPassword()=" << mfc->getPassword() << endl;
  cout << "mfc->getValve()=" << mfc->getValve() << endl;
  cout << "mfc->getOutput()=" << mfc->getOutput() << endl;
  cout << "mfc->getUnits()=" << mfc->getUnits() << endl;
  cout << "mfc->getMax()=" << mfc->getMax() << endl;
  cout << "mfc->getGas()=" << mfc->getGas() << endl;

  /*
  for(Int_t i=0; i<10; ++i){
    Double_t pause=0.+1./10.*(i+1.);
    mfc->setPause(pause);

    Int_t thisPauseSuccess=0;
    Int_t thisPauseTotal=0;
    for(Int_t j=0; j<10; ++j){
      
      ++thisPauseTotal;
      Double_t value=mfc->getGas();
      if(j==0) cout << "value=" << value << endl;
      if( value>0 ) ++thisPauseSuccess;
    }
    Double_t successRate=thisPauseSuccess/thisPauseTotal;
    cout << "pause=" << pause << " " << successRate*100. << "%" 
	 << endl;
  }
  */

  // make sure we close the connection to the serial device
  mfc->Delete();
}

