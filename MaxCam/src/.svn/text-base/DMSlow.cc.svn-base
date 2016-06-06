#include <iostream>
using std::cout;
using std::endl;

#include "TROOT.h"
#include "TSystem.h" 
#include "MaxCamChannel.hh"
#include "TThread.h" 
#include "TMutex.h" 
#include <cmath>

MaxCamChannel * pressure; 
MaxCamChannel * wire_i; 
MaxCamChannel * mesh_i; 
MaxCamChannel * mesh_hv; 
MaxCamChannel * wire_hv;
MaxCamChannel * temp0; 
MaxCamChannel * humidity; 
MaxCamChannel * oxygen; 

TThread * ramp_thr; 
TThread * mesh_thr; 
TThread * wire_thr; 


TMutex * mesh_mutex; 
TMutex * wire_mutex; 


void readPressure() {
  pressure->readFromDB("pressure");
  pressure->readSerialInficon();   
  //cout << pressure->GetName() << ": " << pressure->currentValue   << "  +/- " << pressure->currentValueRMS   <<  " set=" << pressure->setValue << endl;
  static double oldValue=-1;
  if (fabs(pressure->currentValue-oldValue)/pressure->currentValue>1e-3) { // save only changes!!
    double tmp=pressure->currentValue;
    if (oldValue>0) {
      pressure->currentValue=oldValue;
      //pressure->saveValueToDBAndCheck("pressure");
    }
    pressure->currentValue=tmp;
    pressure->saveValueToDBAndCheck("pressure");
    oldValue=pressure->currentValue;
  }
  //pressure->cleanOldValues();
}

void readAnodeCurrent() {
  wire_i->readADC();   
  wire_i->currentValue/=10; 
  wire_i->currentValueRMS/=10; 
  //cout << wire_i->GetName() << ": " << wire_i->currentValue   << "  +/- " << wire_i->currentValueRMS   <<  " set=" << wire_i->setValue << endl;
  static double oldValue=-1;
  if (fabs(wire_i->currentValue-oldValue)>-1*2e-2) { // save only changes!!
    double tmp=wire_i->currentValue;
//     if (oldValue>0) {
//       wire_i->currentValue=oldValue;
//       wire_i->saveValueToDB("wire_i");
//     }
    wire_i->currentValue=tmp;
    wire_i->saveValueToDB("wire_i");
    oldValue=wire_i->currentValue;
  }
  //wire_i->cleanOldValues("wire_i", 100);
}

void readDriftCurrent() {
  mesh_i->readADC();   
  //cout << ch.GetName() << ": " << ch.currentValue   << "  +/- " << ch.currentValueRMS   <<  " set=" << ch.setValue << endl;
  static double oldValue=-1;
  if (fabs(mesh_i->currentValue-oldValue)>-1e-5) { // save only changes!!
    double tmp=mesh_i->currentValue;
  //  if (oldValue>0) {
  //    mesh_i->currentValue=oldValue;
  //    mesh_i->saveValueToDB("mesh_i");
  //  }
    mesh_i->currentValue=tmp;
    mesh_i->saveValueToDB("mesh_i");
    oldValue=mesh_i->currentValue;
  }
  //mesh_i->cleanOldValues("mesh_i", 100);
}
void setMeshAndSave(double val) {
  mesh_hv->readFromDB("mesh_hv");
  mesh_hv->setValue = val;
  mesh_hv->saveValueToDBAndCheck("mesh_hv");
}

void setWireAndSave(double val) {
  wire_hv->readFromDB("wire_hv");
  wire_hv->setValue = val;
  wire_hv->saveValueToDBAndCheck("wire_hv");
}


void readHumidity()
{
   double b = 0.874086;
   double m = 0.030895068;
  
   humidity->readFromDB("humidity"); 
   double old = humidity->currentValue; 
   humidity->readADC(); 
   humidity->currentValue  = (humidity->currentValue - b)/m; 
   humidity->currentValueRMS  = (humidity->currentValueRMS) /m; 
   if (fabs(humidity->currentValue - old) > 0.1)
   {
      humidity->saveValueToDB("humidity"); 
   }
}


void readOxygen()
{
  oxygen->readFromDB("oxygen"); 
  double old = oxygen->currentValue;
  oxygen->readADC(); 
  oxygen->currentValue*=15.; 
  oxygen->currentValueRMS*=15.; 

  if (fabs(oxygen->currentValueRMS - old) > 0.1)
  {
      oxygen->saveValueToDB("oxygen"); 
  }
}

void readTemperature() {
  double deltaT=-1*0.25;

  static double oldT0=-1;
  temp0->readFromDB("temp0");
  temp0->readADC();
  temp0->currentValue*=100; // V->C
  temp0->currentValueRMS*=100; // V->C
  if (fabs(temp0->currentValue-oldT0)>deltaT)  { temp0->saveValueToDB("temp0"); oldT0=temp0->currentValue; }
  
  //temp0->cleanOldValues(); 
  /*
  static double oldT1=-1;
  MaxCamChannel temp1("temp1","", 5,-1,"/etc/dbaccess.txt");
  temp1.readFromDB("temp1");
  temp1.readADC();
  temp1.currentValue*=100; // V->C
  temp1.currentValueRMS*=100; // V->C
  if (fabs(temp1.currentValue-oldT1)>deltaT)  { temp1.saveValueToDB("temp1"); oldT1=temp1.currentValue; }

  static double oldT2=-1;
  MaxCamChannel temp2("temp2","", 6,-1,"/etc/dbaccess.txt");
  temp2.readFromDB("temp2");
  temp2.readADC();
  temp2.currentValue*=100; // V->C
  temp2.currentValueRMS*=100; // V->C
  if (fabs(temp2.currentValue-oldT2)>deltaT)  { temp2.saveValueToDB("temp2"); oldT2=temp2.currentValue; }
  */

}


//Insures only one instance of the thread will run. If it is running already, don't try again. 
void singletonRun(TThread * thr)
{
  int state = thr->GetState(); 
  if (state == 6 || state == 1)
  {
    if (state == 6)
    {
      thr->Kill(); 
      thr->Join(); 
    }
    thr->Run(); 
  }
}


void rampFromDB()
{
  singletonRun(ramp_thr); 
}

void * rampFromDBImpl( void * args) {
  double deltaHV=0.005; // 5V

  
  // change vol if set values different
  // or, update current value in database

  MaxCamChannel hvstatus("hvstatus","",-1,-1,"/etc/dbaccess.txt");
  hvstatus.readFromDB("hvstatus");
  bool hvON = hvstatus.setValue>0 ? true : false;

  cout << "HV STATUS: " << hvON << endl; 

  if (!hvON) {
    mesh_mutex->Lock();
    mesh_hv->ramp(0.0);
    double last_val = mesh_hv->currentValue;
    mesh_hv->readValueFromDB("mesh_hv",2); 
    mesh_hv->currentValue = last_val; 
    mesh_hv->saveValueToDB("mesh_hv"); 
    mesh_hv->print();
    mesh_mutex->UnLock();

    wire_mutex->Lock();
    wire_hv->ramp(0.0);
    last_val = wire_hv->currentValue;
    wire_hv->readValueFromDB("wire_hv",2); 
    wire_hv->currentValue = last_val; 
    wire_hv->saveValueToDB("wire_hv"); 
    wire_hv->print();
    wire_mutex->UnLock();
    //mesh_hv->cleanOldValues(); 
    //wire_hv->cleanOldValues(); 
    return 0;
  }


  mesh_mutex->Lock();
  static double mesh_old_set=-1;
  double mesh_set = mesh_hv->readValueFromDB("mesh_hv",2);
  double mesh_old_val = mesh_hv->readValueFromDB("mesh_hv",0);
  mesh_hv->readADC(1000);
  mesh_hv->print();
  if (fabs(mesh_old_set-mesh_hv->currentValue)>1e-4) { 
    mesh_hv->rampFromDB("mesh_hv"); 
    mesh_old_set=mesh_set; 
  } else if (fabs(mesh_old_val-mesh_hv->currentValue)>deltaHV) { 
    mesh_hv->saveValueToDB("mesh_hv"); 
    mesh_old_val=mesh_hv->currentValue; 
  }
  mesh_mutex->UnLock();

  gSystem->Sleep(500);


  wire_mutex->Lock();
  static double wire_old_set=-1;
  double wire_set = wire_hv->readValueFromDB("wire_hv",2);
  double wire_old_val = wire_hv->readValueFromDB("wire_hv",0);
  wire_hv->readADC(1000);
  wire_hv->print();
  if (fabs(wire_old_set-wire_hv->currentValue)>1e-4)      { 
    wire_hv->rampFromDB("wire_hv"); 
    wire_old_set=wire_set; 
  }
  else if (fabs(wire_old_val-wire_hv->currentValue)>deltaHV) { 
    wire_hv->saveValueToDB("wire_hv"); 
    wire_old_val=wire_hv->currentValue; 
  }
  wire_mutex->UnLock();
  gSystem->Sleep(500);
  //mesh_hv->cleanOldValues(); 
  //wire_hv->cleanOldValues(); 
  return 0; 

}

void *readMeshHV_thr(void*args)
{
  mesh_mutex->Lock();
  mesh_hv->readADC(); 
  mesh_hv->print(); 
  mesh_hv->saveValueToDB("mesh_hv",false); 
  mesh_mutex->UnLock();
  return 0;
}

void *readWireHV_thr(void*args)
{
  wire_mutex->Lock();
  wire_hv->readADC(); 
  wire_hv->print(); 
  wire_hv->saveValueToDB("wire_hv",false); 
  wire_mutex->UnLock();
  return 0;
}

void readVoltages()
{
  singletonRun(mesh_thr);   
  singletonRun(wire_thr);   
}
 

int main(int argc, char *argv[]) {


  pressure = new MaxCamChannel("pressure"); 
  wire_i = new MaxCamChannel("wire_i"); 
  mesh_i = new MaxCamChannel("mesh_i"); 
  mesh_hv = new MaxCamChannel("mesh_hv"); 
  wire_hv = new MaxCamChannel("wire_hv"); 
  temp0 = new MaxCamChannel("temp0"); 
  humidity = new MaxCamChannel("humidity"); 
  oxygen = new MaxCamChannel("oxygen"); 

  ramp_thr = new TThread("ramp_thr",rampFromDBImpl,NULL); 
  mesh_thr = new TThread("mesh_thr",readMeshHV_thr,NULL); 
  wire_thr = new TThread("wire_thr",readWireHV_thr,NULL); 


  TMutex * mut_db = new TMutex(); 
  TMutex * mut_adc = new TMutex(); 

  wire_hv->synchronizeComedi(mut_adc); 
  wire_i->synchronizeComedi(mut_adc); 
  mesh_hv->synchronizeComedi(mut_adc); 
  mesh_i->synchronizeComedi(mut_adc); 
  temp0->synchronizeComedi(mut_adc); 
  oxygen->synchronizeComedi(mut_adc); 
  humidity->synchronizeComedi(mut_adc); 

  wire_hv->synchronizeDB(mut_db); 
  pressure->synchronizeDB(mut_db); 
  wire_i->synchronizeDB(mut_db); 
  mesh_hv->synchronizeDB(mut_db); 
  mesh_i->synchronizeDB(mut_db); 
  temp0->synchronizeDB(mut_db); 
  humidity->synchronizeDB(mut_db); 
  oxygen->synchronizeDB(mut_db); 
  
  mesh_mutex = new TMutex(); 
  wire_mutex = new TMutex(); 

//  mesh_hv->calibrateOffset(); 
//  wire_hv->calibrateOffset(); 

  int nsleep=500;
  int nskip=1; // read voltage every nsleep*mskip msec
  if (argc>1) {
	nsleep=atoi( argv[1] );
  	cout << argc << "   " << (argv[1]) << endl;
  }
  cout << "Sleep set to " << nsleep << " ms" << endl;

  int nread=0;
  while (1) {
    readPressure();
    if (nread++ % nskip == 0) {
        
//        readVoltages(); 
        readAnodeCurrent();
        readDriftCurrent();
        rampFromDB();
        readTemperature();
        readHumidity(); 
        readOxygen(); 
    }	
    gSystem->Sleep(nsleep);
  }

  return 0;
}
