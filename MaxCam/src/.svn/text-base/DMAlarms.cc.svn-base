#include <iostream>
using std::cout;
using std::endl;

#include "TROOT.h"
#include "TSystem.h" 
#include "TSocket.h"
#include "MaxCamChannel.hh"
#include <cmath>
#include <unistd.h>


MaxCamChannel * pressure; 
MaxCamChannel * wire_i; 
MaxCamChannel * wire_hv;
MaxCamChannel * mesh_hv;
MaxCamChannel * ups;
MaxCamChannel * hvstatus;

double pressureAlert = 2.0;
double pressureAlarm = 5.0;

double anodeIAlert = 0.02;
double anodeIAlarm = 0.03;

double anodeVAlert = 0.02;
double anodeVAlarm = 0.05;

double driftVAlert = 0.05;
double driftVAlarm = 0.1;

//Return values of functions: 1=OK, 0=Alert, -1=Alarm

double lastAlertPressure =0;
double lastAlertAnodeI=0;
double lastAlertAnodeV=0;
double lastAlertSlowControl=0;
double lastAlertDAQ=0;
double lastAlertPower=0;
double lastAlertNetwork=0;
double lastAlertDriftV=0;

/** Sends an alarm to the system \param messagetext text of the alarm to be sent*/
void sendAlarm(TString messagetext)
{
  TString hnmessage = "\""+messagetext+"\"";
  cout << messagetext << endl;
  system("./send_alert \"!!! 10L ALARM !!!\" "+hnmessage);
  system("mysql DM_SLOWCONTROL -u dmatter -pseedark -e \"insert into alerts (alarm,content,timestamp) VALUES(1,'" + hnmessage + "',NOW())\""); 


}

/** Sends an alert to the system \param messagetext text of the alert to be sent*/
void sendAlert(TString messagetext)
{
  TString hnmessage = "\""+messagetext+"\"";
  cout << messagetext << endl;
  system("./send_alert \"10L Alert\" "+hnmessage);
  system("mysql DM_SLOWCONTROL -u dmatter -pseedark -e \"insert into alerts (alarm,content,timestamp) VALUES(0,'" + hnmessage + "',NOW())\""); 
}

int checkPressure()
{
  pressure->readFromDB("pressure");
  double setP = pressure->setValue;
  double currentP = pressure->currentValue;

  cout << "setP: " << setP << " currentP: " << currentP << endl;

  bool isBusy = pressure->isBusy();


  if(fabs(currentP-setP) > pressureAlarm && !isBusy)
    {
      TString text = "Pressure has triggered an alarm! Current pressure is: ";
      text+=currentP;
      text+=". Set Pressure is: ";
      text+=setP;
      text+=". Beginning emergency shut down!";
      sendAlarm(text);
      
      return -1;
    }
  else if(fabs(currentP-setP) > pressureAlert && !isBusy && time(NULL)-lastAlertPressure > 1800)
    {
      TString text = "Pressure has triggered an alert! Current pressure is: ";
      text+=currentP;
      text+=". Set Pressure is: ";
      text+=setP;
      text+=".";
      sendAlert(text);
      lastAlertPressure=time(NULL);
	
      return 0;
    }
  else
    {
      cout << "Pressure OK!" << endl;
      return 1;
    }
}

int checkAnode()
{
  wire_i->readFromDB("wire_i");
  wire_hv->readFromDB("wire_hv");
  mesh_hv->readFromDB("mesh_hv");
  double currentI =  wire_i->getAverage(10);
  double currentV = wire_hv->getAverage(10);
  double setV = wire_hv->setValue;
  
  hvstatus->readFromDB("hvstatus");
  bool hvON = hvstatus->setValue>0 ? true : false;
  
  bool isBusyWire = wire_hv->isBusy();
  bool isBusyMesh = mesh_hv->isBusy();
  bool isBusy = (isBusyWire || isBusyMesh);
  
  if(currentI > anodeIAlarm && !isBusy)
    {
      TString text = "Anode current has triggered an alarm! Anode current is ";
      text+=currentI;
      text+=". Beginning emergency shut down!";
      sendAlarm(text);
      return -1;
    }
  else if(currentI > anodeIAlert && !isBusy && time(NULL)-lastAlertAnodeI)
    {
      TString text="Anode current has triggered an alert! Anode current is ";
      text+=currentI;
      text+=".";
      sendAlert(text);
      lastAlertAnodeI=time(NULL);
      return 0;
    }
//   else if (fabs(setV-currentV)> anodeVAlarm && !isBusy && hvON)
//     {
//       TString text = "Anode voltage has triggered an alarm! Anode voltage is: ";
//       text+=currentV;
//       text+=". The set voltage is: ";
//       text+=setV;
//       text+=". Beginning emergency shutdown!";
//       sendAlarm(text);
//       return -1;
//     }
//   else if( fabs(setV-currentV)> anodeVAlert && !isBusy && hvON && time(NULL)-lastAlertAnodeV >1800)
//     {
//       TString text = "Anode voltage has triggered an alarm! Anode voltage is: ";
//       text+=currentV;
//       text+=". The set voltage is: ";
//       text+=setV;
//       text+=".";
//       sendAlert(text);
//       lastAlertAnodeV=time(NULL);
//       return 0;
//     }
  else
    {
      cout << "Anode OK!" << endl;
      return 1;
    }
}

int checkDrift()
{
  mesh_hv->readFromDB("mesh_hv");
  double currentV = mesh_hv->getAverage(10);
  double setV = mesh_hv->setValue;
  
  hvstatus->readFromDB("hvstatus");
  bool hvON = hvstatus->setValue>0 ? true : false;
  
  bool isBusyWire = wire_hv->isBusy();
  bool isBusyMesh = mesh_hv->isBusy();
  bool isBusy = (isBusyWire || isBusyMesh);
  
  if (fabs(setV-currentV)> driftVAlarm && !isBusy && hvON)
    {
      TString text = "Drift voltage has triggered an alarm! Drift voltage is: ";
      text+=currentV;
      text+=". The set voltage is: ";
      text+=setV;
      text+=". Beginning emergency shutdown!";
      sendAlarm(text);
      return -1;
    }
  else if( fabs(setV-currentV)> driftVAlert && !isBusy && hvON && time(NULL)-lastAlertDriftV >1800)
    {
      TString text = "Drift voltage has triggered an alarm! Drift voltage is: ";
      text+=currentV;
      text+=". The set voltage is: ";
      text+=setV;
      text+=".";
      sendAlert(text);
      lastAlertDriftV=time(NULL);
      return 0;
    }
  else
    {
      cout << "Drift OK!" << endl;
      return 1;
    }
}

int checkSlowControl()
{
  TString checkslowlines = gSystem->GetFromPipe("ps aux | awk /DMSlow$/ | wc -l");
  TString checkDB = gSystem->GetFromPipe("mysql DM_SLOWCONTROL -u dmatter -pseedark -e \"select dmslow from busy\" | awk /[0-1]/");

  if(checkslowlines>=checkDB)
    {
      cout << "Slow control OK!" << endl;
      return 1;
    }
  else if(time(NULL)-lastAlertSlowControl >1800)
    {
      sendAlert("Check Slow Control Process!");
      lastAlertSlowControl=time(NULL);
      return 0;
    }
  
  
  return 0;
  
}

int checkDAQ()
{
 
  TString daqlines = gSystem->GetFromPipe("ssh dmatter@gimli.lns.mit.edu 'ps ax | awk \"\\/usr\\/bin\\/python \\.\\/run10L\\/\" | wc -l'"); 
  TString daqdb = gSystem->GetFromPipe("mysql DM_SLOWCONTROL -u dmatter -pseedark -e \"select daq from busy\" | awk /[0-1]/"); 
  TString ccdinitlines = gSystem->GetFromPipe("ssh dmatter@gimli.lns.mit.edu 'ps ax | awk \"\\/bin\\/bash.*ccdinit\\/\" | wc -l'"); 
  TString daqsublines  = gSystem->GetFromPipe("ssh dmatter@gimli.lns.mit.edu 'ps ax | awk \"\\/root .* run10L.cxx\([0-9]/\"  | wc -l '"); 
  TString daqsubdb = gSystem->GetFromPipe("mysql DM_SLOWCONTROL -u dmatter -pseedark -e \"select daq_sub from busy\" | awk /[0-1]/"); 

  if (daqlines.Atoi() + ccdinitlines.Atoi() >= daqdb.Atoi() && daqsublines.Atoi() >= daqsubdb.Atoi() )
  {
      cout << "DAQ  OK!" << endl;
      return 1;
  }
  else if(time(NULL)-lastAlertDAQ >1800)
    {
      sendAlert("Check DAQ!");
      lastAlertDAQ=time(NULL);
      return 0;
    }


  return 1; 

}

int checkPower()
{
  int treebeard = atoi(gSystem->GetFromPipe("./DM_UPS treebeard.lns.mit.edu > /dev/null; echo $?")); 
  eout << " Treebeard returns: " << treebeard << endl; 
  int boromir = atoi(gSystem->GetFromPipe("./DM_UPS boromir.lns.mit.edu > /dev/null; echo $?")); 
  cout << " Boromir returns: " << boromir << endl; 

  if (!treebeard && !boromir)
  {
    cout << "Power on!" << endl; 
    return 1; 
  }
  if (treebeard == 10 && boromir == 0)
  {
    cout << "Treebeard Busy! " << std::endl; 
    return 1; 
  }
  if (treebeard == 0 && boromir == 10)
  {
    cout << "Boromir Busy! " << std::endl; 
    return 1; 
  }
  
  if (treebeard == 10 && boromir == 10)
  {
    cout << " Both Busy!" << endl; 
    return 1; 
  }

  if (treebeard == 1 || treebeard == 2 || treebeard == 3 || boromir == 1 || boromir == 2 || boromir ==3)
  {
    sendAlert( "Power off!  Beginning Emergency Shutdown!"); 
    return -1; 
  }

  sendAlert("Error in power system"); 
  return 0; 

}


int checkNetwork()
{
  TSocket* gandalf = new TSocket("gandalf.lns.mit.edu",22);
  TSocket* google = new TSocket("google.com",80);
  
  if(!gandalf->IsValid())
    {
      TString text= "gandalf could not be reached!";
      if(time(NULL)-lastAlertNetwork>1800)
      {
        sendAlert(text);
        lastAlertNetwork=time(NULL);
      }
      gandalf->Close();
      delete gandalf;
      return -1;
    }
  
   if(!google->IsValid())
    {
      TString text= "google could not be reached!";
      if(time(NULL)-lastAlertNetwork>1800)
      {
        sendAlert(text);
        lastAlertNetwork=time(NULL);
      }
      google->Close();
      delete google;
      return -1;
    }
  
  cout << "Network up!" << endl;
  google->Close();
  delete google;
  gandalf->Close();
  delete gandalf;
  return 1;

}

int main(int argc, char *argv[]) {

  pressure = new MaxCamChannel("pressure"); 
  wire_i = new MaxCamChannel("wire_i"); 
  wire_hv = new MaxCamChannel("wire_hv"); 
  mesh_hv = new MaxCamChannel("mesh_hv");
  hvstatus = new MaxCamChannel("hvstatus");

  int nsleep=0;
  int nskip=1; // read voltage every nsleep*mskip msec
  if (argc>1) {
    nsleep=atoi( argv[1] );
    cout << argc << "   " << (argv[1]) << endl;
  }
  cout << "Sleep set to " << nsleep << " ms" << endl;

  int nread=0;
  while (1) {
    if(checkPressure()<0)
      break;
    if(checkAnode()<0)
      break;
//     if(checkDrift()<0)
//       break;
    if(checkSlowControl()<0)
      break;
//    if(checkDAQ()<0)
//      break;
    if(checkPower()<0)
      break;
    if(checkNetwork()<0)
      break;
    gSystem->Sleep(nsleep);
    
  }	
  
  cout << "Emergency Condition reached!" << endl;
  system("./EmergencyShutDown");
  
  return 0;
}
