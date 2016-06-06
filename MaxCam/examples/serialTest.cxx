#include <string>
#include <iostream>
#include "../MaxCamSerial.hh"
#include "../MaxCamChannel.hh"
using namespace std;

void testConvectron(char *port="/dev/ttyS3") { 
  //SerialHandle rs232 = MaxCamSerial::Open(port);
  //MaxCamSerial::Init(rs232, 19200, 8, 'n', 1);
  ////MaxCamSerial::Init(rs232);
  //MaxCamSerial::Send(rs232, "RD\r", 0.1);
  //int readlength = 16;
  //string readstr = MaxCamSerial::ReadPort(rs232, readlength);
  //MaxCamSerial::Close(rs232);
  //cout << "readstr = [" << readstr << "]" << endl << endl;

  MaxCamChannel *ser = new MaxCamChannel("convectron", "convectron");
  double pressure = ser->readConvectron(port);
  cout << "pressure = " << pressure << " torr" << endl;
}

void testInficon(char *port_name="/dev/ttyS1", int channel=1) {

  MaxCamChannel *ser = new MaxCamChannel("inficon", "inficon");
  double pressure = ser->readInficonController(port_name, channel);
  cout << "pressure = " << pressure << " torr" << endl;
}
