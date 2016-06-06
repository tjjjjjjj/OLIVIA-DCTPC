#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>  // sleep()
//#include "../MaxCamSerial.hh"
#include "../MaxCamChannel.hh"
#include "TTimeStamp.h"
using namespace std;

// void testConvectron(char *port="/dev/ttyS3") { 
//   //SerialHandle rs232 = MaxCamSerial::Open(port);
//   //MaxCamSerial::Init(rs232, 19200, 8, 'n', 1);
//   ////MaxCamSerial::Init(rs232);
//   //MaxCamSerial::Send(rs232, "RD\r", 0.1);
//   //int readlength = 16;
//   //string readstr = MaxCamSerial::ReadPort(rs232, readlength);
//   //MaxCamSerial::Close(rs232);
//   //cout << "readstr = [" << readstr << "]" << endl << endl;

//   MaxCamChannel *ser = new MaxCamChannel("convectron", "convectron");
//   double pressure = ser->readConvectron(port);
//   cout << "pressure = " << pressure << " torr" << endl;
// }

// void testInficon(char *port_name="/dev/ttyS1", int channel=1) {

//   MaxCamChannel *ser = new MaxCamChannel("inficon", "inficon");
//   double pressure = ser->readInficonController(port_name, channel);
//   cout << "pressure = " << pressure << " torr" << endl;
//}

void PressuretoFile(int nsamples=1, float samplrate=1, char* inficport="/dev/ttyS1", char *convecport="/dev/ttyS3", int channel=1, char* outputfile="pressurefile.dat") {
  
  ofstream myfile;
  myfile.open(outputfile);
  myfile<<"# time(yyyymmdd-hhmmss) inficon["<<inficport<<"](torr) convectron["<<convecport<<"](torr) \n";
 
  MaxCamChannel *cser = new MaxCamChannel("convectron", "convectron");
  MaxCamChannel *iser = new MaxCamChannel("inficon", "inficon");

  for (int n=0; n<nsamples; n++){
    //cout << "n = " << n << endl;
    double cpressure = cser->readConvectron(convecport);  
    //cout << "convectron done" << endl;
    double ipressure = iser->readInficonController(inficport, channel);
    //cout << "inficon done" << endl;
      
    TTimeStamp time;

    myfile<<time.GetDate()<<"-"<<time.GetTime()<<" "<<ipressure<<" "<<cpressure<<endl;
    //cout << "file io done" << endl;
      
    if (samplrate > 0.0) {
      //cout << "sleep for " << samplrate << endl;
      sleep(samplrate);
    }
    //cout << "sleep done" << endl;
  }

  cout << "Closing " << outputfile << endl;
  myfile.close();
}
