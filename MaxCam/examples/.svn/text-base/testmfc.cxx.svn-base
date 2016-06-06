#include <iostream>
#include "../MaxCamSerial.hh"
#include "../MaxCamMFC.hh"
#include <climits>
#include <sstream>
#include <iomanip>  // setprecision(), fixed
using namespace std;

volatile int char_min = CHAR_MIN;

void testread() {
  // open up connection to MFC
  SerialHandle mfc = MaxCamSerial::Open("/dev/ttyS2");
  MaxCamSerial::Init(mfc, 9600, 8, 'n', 1);
  int readlength = 255;
  bool trimWhitespace = true;
  string readstr = MaxCamSerial::ReadPort(mfc, readlength, trimWhitespace);
  cout << "readstr = " << readstr << endl;

  cout << endl << " ... " << endl;

  // get rid of anything after the last carriage return '\r'
  int endpos = readstr.find_last_of('\r');
  cout << "endpos = " << endpos << endl;
  readstr = readstr.substr(0, endpos);
  cout << "readstr = " << readstr << endl;
}

void test() {

  // open up connection to MFC
  SerialHandle mfc = MaxCamSerial::Open("/dev/ttyS2");
  MaxCamSerial::Init(mfc, 9600, 8, 'n', 1);
  int readlength = 255;
  bool trimWhitespace = true;
  //MaxCamSerial::Flush(mfc);

  char sendstring[8];
  sendstring[0] = '?';
  sendstring[1] = 'S';
  sendstring[2] = 'i';
  sendstring[3] = 'n';
  sendstring[4] = 'v';
  sendstring[5] = '\xA5';
  sendstring[6] = '\x72';
  sendstring[7] = '\0';

  unsigned int crc = MaxCamMFC::CalcCRC("?Sinv");
  cout << hex << "CRC = " << crc << endl;

  //char sendstring[13];
  //sendstring[0]  = 'S';
  //sendstring[1]  = 'i';
  //sendstring[2]  = 'n';
  //sendstring[3]  = 'v';
  //sendstring[4]  = '2';
  //sendstring[5]  = '.';
  //sendstring[6]  = '0';
  //sendstring[7]  = '0';
  //sendstring[8]  = '0';
  //sendstring[9]  = '\x8F';
  //sendstring[10] = '\x55';
  //sendstring[11] = '\r';
  //sendstring[12] = '\0';

  MaxCamSerial::Send(mfc, sendstring, 0.1);
  string readstr = MaxCamSerial::ReadPort(mfc, readlength, trimWhitespace);
  cout << "readstr = " << readstr << endl;

}

void makeMFCSendString(string macro, float val) {

  char letter;
  stringstream cmdstream;
  //cmdstream << "Sinv2.000";
  //cmdstream << "?Val";

  //cmdstream << "Sinv";
  cmdstream << macro;
  //cmdstream.setf(ios::fixed, ios::floatfield);
  //cmdstream.precision(ndecimalPlaces);  
  cmdstream << std::fixed;
  int ndecimalPlaces = 3;
  cmdstream << std::setprecision(ndecimalPlaces);
  cmdstream << val;

  //cmdstream << "48";  // hex code for "H"
  cout   << cmdstream.str() << endl;
  cout << "cmdstream.str().length() = " << cmdstream.str().length() << endl;

  int nbytesAppend = 3;  // 2 CRC bytes and a carriage return
  int stringLength = cmdstream.str().length();
  const int chararrLength = stringLength + nbytesAppend;
  cout << "chararrLength = " << chararrLength << endl;

  // copy from the string stream to the chararr
  // is there a better way?
  char chararr[chararrLength];
  for (int ii=0; ii<stringLength; ii++) {
    cmdstream >> letter;
    cout << "letter = " << letter << endl;
    chararr[ii] = letter;
  }

  // Compute and append the CRC bits
  char cmdForCrc[chararrLength];
  strcpy(cmdForCrc, cmdstream.str().c_str());
  unsigned int crc = MaxCamMFC::CalcCRC(cmdForCrc);
  cout << "CRC = " << hex << crc << endl;

  int crcHigh = (crc >> 8) & 0xFF;
  int crcLow  = crc & 0xFF;
  cout << "crcHigh = " << hex << crcHigh << endl;
  cout << "crcLow  = " << hex << crcLow  << endl;
  // bits for Sinv2.000
  //int crcHigh = 0x8F;
  //int crcLow  = 0x55;
  // bits for ?Val
  //int crcHigh = 0x05;
  //int crcLow  = 0x98;
  chararr[stringLength]   = crcHigh;
  chararr[stringLength+1] = crcLow;
  
  // add the carriage return and null-terminate
  chararr[stringLength+2] = '\r';
  chararr[chararrLength] = '\0';  
  cout << "chararr = " << chararr << endl;

  for (int ii=0; ii<chararrLength; ii++) {
    cout << "char = " << chararr[ii] << " = " << hex << (int)chararr[ii] << dec << endl;
  }

}

void test2() {

  unsigned int crc = 0x8F55;
  const int CMDLENGTH = 13;

  char          cmd [CMDLENGTH];
  unsigned char cmdu[CMDLENGTH];

  int          crcTop  = 0x8F;
  unsigned int crcTopU = 0x8F;

  cmd[0]  = (unsigned int)crcTop;
  cmdu[0] = crcTopU;
  cout << hex << showbase;
    
  cout << "crcTop  = " << crcTop << endl;
  cout << "crcTopU = " << crcTopU << endl;
  cout << "(int)cmd[0]           = " << (int)cmd[0] << endl;
  cout << "(unsigned int)cmd[0]  = " << (unsigned int)cmd[0] << endl;
  cout << "     cmd[0]           = " << cmd[0] << endl;
  cout << "(unsigned int)cmdu[0] = " << (unsigned int)cmdu[0] << endl;
  cout << "              cmdu[0] = " << cmdu[0] << endl;

  if (char_min < 0)
    cout << "char is signed";
  else if (char_min == 0)
    cout << "char is unsigned";
  else
    cout << "unrecognized char";
  cout << endl;
}

void mfc() {
  
  const int CMDLENGTH=13;
  char cmd[CMDLENGTH];
  strcpy(cmd, "Sinv2.000");
  unsigned int crc = MaxCamMFC::CalcCRC(cmd);
  cout << hex << crc << endl;

  // get the two low bytes and the high two bytes
  // e.g. if crc = 0x8F55 then crcDown = 0x55 and crcUp = 0x8F
  unsigned int crcDown = crc & 0x00FF;
  //unsigned int crcUp   = (crc>>8);  // shift right by two slots
  unsigned int crcUp   = ( (crc & 0xFFFF) & 0xFF00)>>8;  // shift right by one byte

  //unsigned int crcUp   = (crc>>8) & 0x00FF;  // shift right by two slots
  cout << "hi: [" << hex << crcUp   << "] = " << (char)crcUp   << "]" << endl;
  cout << "lo: [" << hex << crcDown << "] = " << (char)crcDown << "]" << endl;

  //cmd[9]  = (char)crcUp;
  cmd[9]  = crcUp;
  cmd[10] = crcDown;
  cmd[11] = '\r';
  cmd[12] = '\0';
  cout << showbase;
  cout << hex << "cmd[9]  = [" << cmd[9] << "]" << endl;
  cout << hex << "cmd[9]  = [" << (int)cmd[9] << "]" << endl;
  cout << hex << "cmd[10] = [" << cmd[10] << "]" << endl;
  cout << hex << "cmd[10] = [" << (int)cmd[10] << "]" << endl;
  cout << dec;
  for (int ii=0; ii<CMDLENGTH; ii++) {
    //cout << hex << (int)cmd[ii] << " ";
    cout << ii << " = " << hex << cmd[ii] << " = " << (int)cmd[ii] << dec << " " << endl;;
  }
  cout << endl;





  //SerialHandle mfc = MaxCamSerial::Open("/dev/ttyS2");
  //MaxCamSerial::Init(mfc, 9600, 8, 'n', 1);
  //
  //int readlength = 24;
  ////string readstr = MaxCamSerial::ReadPort(mfc, readlength);
  //char readbuf[64];
  //char *readptr;
  //for (int ii=0; ii<64; ii++) {
  //  readbuf[ii] = 'x';
  //}
  //readptr = readbuf;
  //readptr = MaxCamSerial::ReadPort(mfc, readbuf, readlength);
  //MaxCamSerial::Close(mfc);
  //
  ////cout << "readstr = [" << readstr << "]" << endl;
  //cout << "readbuf = [" << readbuf << "]" << endl;

}
