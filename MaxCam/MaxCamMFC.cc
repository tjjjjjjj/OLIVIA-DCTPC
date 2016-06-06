

#include "MaxCamMFC.hh"

#include <string> 
#include <cstdlib>
#include <termios.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <unistd.h>
//#include <fstream>
//#include <vector>
//#include <iomanip>
//#include <termios.h>
//#include <fcntl.h>
//#include <unistd.h>

//#include <math.h>

//#ifdef DM_DAQ
//#endif

//namespace DBAccess {
//  char server[256];
//  char user[256];
//  char pass[256];
//  char database[256];
//  void set(char *fname) {
//    ifstream infile(fname);
//    infile >> server >> user >> pass >> database;
//  }
//};

// Cyclic Redundancy Check codes for 
// the Sierra Instruments command set.
// command crc_code
// ?Ver    0x3aa3
// ?Ffs    0x3cb2
// ?Srn    0x980a
// ?Sinv   0xa572
// ?Sini   0x46ac
// ?Xwrd   0x9762
// ?Vlvi   0x9bc3
// ?Outi   0xc883
// ?Unti   0x081d
// ?Fscl   0x4315
// ?Gasi   0x4b74

// Vlvi2   0x1250

ClassImp(MaxCamMFC)

  MaxCamMFC::MaxCamMFC(const char *name, const char *title, char *port, string model, char *dbaccess) :
  TNamed(name, title),
  currentFlowRate(0), 
  port(port)
{
  cout << GetName() << endl;
  if(debug)  cout << "opening serial port" << endl;

  // the default way of opening the serial port 
  //  rs232 = MaxCamSerial::Open(port);
  // with all options

  // nope
  rs232 = open(port, O_RDWR | O_NOCTTY | O_NDELAY);

  // noped
  //  rs232 = open(port, O_RDWR | O_NDELAY);
  
  // nope
  //  rs232 = open(port, O_RDWR | O_NOCTTY );

  // nope
  //  rs232 = open(port, O_NOCTTY | O_NDELAY);

  // nope
  //  rs232 = open(port, O_RDWR );

  // nope
  //  rs232 = open(port, O_NOCTTY );

  // nope
  //  rs232 = open(port, O_NDELAY);

  MaxCamSerial::Init(rs232, 9600, 8, 'n', 1);
  if(debug)  cout << "rs232=" << rs232 << endl;

  // set debug flag off by default
  debug=false;

  // pause this much on send commands
  pause=1.;
  
  //cout << "MaxCamMFC ctor" << endl;
  //if (dbparam) DBAccess::set(dbparam);
}

MaxCamMFC::MaxCamMFC(const MaxCamMFC &other) :
  TNamed(other),
  currentFlowRate(other.currentFlowRate),
  baudRate(other.baudRate),
  port(port)
{}

void 
MaxCamMFC::setPause(float p) {
  pause=p;
}

void 
MaxCamMFC::setValve(int ii) {
  // VlviValueCCR
  // Value can be an integer between 1 and 3
  // 1=Automatic, controls the set point value
  // 2=Closed
  // 3=Purge

  char cmd[32];
  sprintf(cmd,"Vlvi%d",ii);
  // append CRC to end of string
  Int_t nbytes=sprintf(cmd,"%s%s",cmd,CalcCRC(cmd).c_str());
  //  cout << "cmd=" << cmd << endl;
  
  MaxCamSerial::Flush(rs232);
  MySerialWrite(rs232,cmd,nbytes);
}

void 
MaxCamMFC::setUnits(int ii) {
  // UntiValueCCR
  // Value can be an integer between 1 and 3
  // 1=scc/s
  // 2=scc/m
  // 3=scc/H
  // 4=Ncc/s
  // 5=Ncc/m
  // 6=Ncc/H
  // 7=SCF/s
  // 8=SCF/m
  // 9=SCF/
  // 10=NM^3/s
  // 11=NM^3/M
  // 12=NM^3/H
  // 13=SM^3/s
  // 14=SM^3/M
  // 15=SM^3/H
  // 16=sl/s
  // 17=sl/M
  // 18=sl/H
  // 19=NL/s
  // 20=NL/M
  // 21=NL/H
  // 22=g/s
  // 23=g/M
  // 24=g/H
  // 25=kg/s
  // 26=kg/M
  // 27=kg/H
  // 28=lb/s
  // 29=lb/M
  // 30=lb/H

  char cmd[32];
  sprintf(cmd,"Unti%d",ii);
  // append CRC to end of string
  Int_t nbytes=sprintf(cmd,"%s%s",cmd,CalcCRC(cmd).c_str());
  //  cout << "cmd=" << cmd << endl;
  
  MaxCamSerial::Flush(rs232);
  MySerialWrite(rs232,cmd, nbytes);
}

void 
MaxCamMFC::setGas(int ii) {
  // GasiValueCCR
  // Value can be an integer between 1 and ?
  // 1=Air
  // 2=Argon
  // 3=Carbon Dioxide
  // 4=Carbon Monoxide
  // 5=Helium
  // 6=90% Argon/10% Carbon Dioxide
  // 7=90% Argon/10% Methane
  // 8=Carbon Tetra-Fluoride
  // 9=Helium 3
  // 10=CS2

  char cmd[32];
  sprintf(cmd,"Gasi%d",ii);
  // append CRC to end of string
  Int_t nbytes=sprintf(cmd,"%s%s",cmd,CalcCRC(cmd).c_str());
  //  cout << "cmd=" << cmd << endl;
  
  MaxCamSerial::Flush(rs232);
  MySerialWrite(rs232,cmd, nbytes);
}

void 
MaxCamMFC::setSource(int ii) {
  // SiniValueCCR
  // Value can be an integer between 1 and 5
  // 1=Module/RS232 as the set point source
  // 2=0-5 V analog set point source
  // 3=0-10 V analog set point source
  // 4=1-5 volt analog set point source
  // 5=4-20 ma analog set point source

  char cmd[32];
  sprintf(cmd,"Sini%d",ii);
  // append CRC to end of string
  Int_t nbytes=sprintf(cmd,"%s%s",cmd,CalcCRC(cmd).c_str());
  //  cout << "cmd=" << cmd << endl;
  
  MaxCamSerial::Flush(rs232);
  MySerialWrite(rs232,cmd, nbytes);
}

void 
MaxCamMFC::setOutput(int ii) {
  // OutiValueCCR
  // Value can be an integer between 1 and 3
  // 1=0-5 volt and 4-20 ma outputs
  // 2=0-10 volt and 4-20 ma outputs
  // 3=1-5 volt and 4-20 ma outputs


  char cmd[32];
  sprintf(cmd,"Outi%d",ii);
  // append CRC to end of string
  Int_t nbytes=sprintf(cmd,"%s%s",cmd,CalcCRC(cmd).c_str());
  //  cout << "cmd=" << cmd << endl;
  
  MaxCamSerial::Flush(rs232);
  MySerialWrite(rs232,cmd, nbytes);
}

string
MaxCamMFC::queryMFC(char* c,int arglength) {
  // ?VerCCR
  char cmd[32],cmdstart[32];
  sprintf(cmd,"%s",c);
  sprintf(cmdstart,"%s",c);
  
  // append CRC to end of command
  Int_t nbytes=sprintf(cmd,"%s%s",cmd,CalcCRC(cmd).c_str());

  //  cout << "cmd=" << cmd << endl;
  
  // flush the input buffer before we send the command
  MaxCamSerial::Flush(rs232);
  MySerialWrite(rs232, cmd, nbytes);

  MaxCamSerial::Wait(pause);

  // if we print the buffer, it clears the buffer

  //  printBuffer();
  //  cout << "cmdstart=" << cmdstart << endl;
  if(debug) cout << "cmdstart=" << cmdstart << endl;
  string answer=GetFeedback(cmdstart,5);
  if(debug)  cout << "answer=" << answer << endl;
  return answer;
}

string
MaxCamMFC::getVersion(void) {

  string version=queryMFC("?Ver",5);
  if(debug)  cout << "version=" << version << endl;
  return version;

  /* the old getVersion(void)...
  // ?VerCCR
  char cmd[32],cmdstart[32];
  sprintf(cmd,"?Ver");
  sprintf(cmdstart,"?Ver");
  
  // append CRC to end of command
  sprintf(cmd,"%s%s",cmd,CalcCRC(cmd).c_str());
  //  cout << "cmd=" << cmd << endl;
  
  // flush the input buffer before we send the command
  MaxCamSerial::Flush(rs232);
  MySerialWrite(rs232, cmd, nbytes);
  
  //  cout << "cmdstart=" << cmdstart << endl;
  string version=GetFeedback(cmdstart,5);
  if(debug)  cout << "version=" << version << endl;
  return version;
  */
}

void 
MaxCamMFC::openValve(void){
  // open the valve
  setValve(1);
}

void 
MaxCamMFC::closeValve(void){
  // open the valve
  setValve(2);
}

void 
MaxCamMFC::purgeValve(void){
  // open the valve
  setValve(3);
}

void 
MaxCamMFC::Fill(float rate, int seconds) {
  // valve should be open already, but make sure
  setFlowPoint(rate);
  sleep(seconds);
  setFlowPoint(0.000);
}

int
MaxCamMFC::getSerial(void){
  string serial=queryMFC("?Srn",7);
  if(debug)  cout << "serial=" << serial << endl;
  return atoi(serial.c_str());
}

float
MaxCamMFC::getFactor(void){
  string factor=queryMFC("?Ffs",7);
  if(debug)  cout << "factor=" << factor << endl;
  return atof(factor.c_str());
}

string 
MaxCamMFC::getFlowPoint(void){
  string flowpoint=queryMFC("?Sinv",5);
  if(debug)  cout << "flowpoint=" << flowpoint << endl;
  return flowpoint;
}

int
MaxCamMFC::getSource(void){
  string source=queryMFC("?Sini",1);
  if(debug)  cout << "source=" << source << endl;
  return atoi(source.c_str());
}

string
MaxCamMFC::getPassword(void){
  string pwd=queryMFC("?Xwrd",4);
  if(debug)  cout << "pwd=" << pwd << endl;
  return pwd;
}

int
MaxCamMFC::getValve(void){
  string valve=queryMFC("?Vlvi",1);
  if(debug)  cout << "valve=" << valve << endl;
  return atoi(valve.c_str());
}

int
MaxCamMFC::getOutput(void){
  string output=queryMFC("?Outi",1);
  if(debug)  cout << "output=" << output << endl;
  return atoi(output.c_str());
}

int
MaxCamMFC::getUnits(void){
  string units=queryMFC("?Unti",1);
  if(debug)  cout << "units=" << units << endl;
  return atoi(units.c_str());
}

float
MaxCamMFC::getMax(void){
  string max=queryMFC("?Fscl",7);
  if(debug)  cout << "max=" << max << endl;
  return atof(max.c_str());
}

int
MaxCamMFC::getGas(void){
  string gas=queryMFC("?Gasi",1);
  if(debug)  cout << "gas=" << gas << endl;
  return atoi(gas.c_str());
}

float
MaxCamMFC::getFlow(void) {
  // SHOULD ADD SOMETHING THAT ACCOUNTS FOR FACT THAT SOMETIMES THERE'S
  // NOTHING IN THE BUFFER YET...
  // FlowValueCCR
  // flush the input buffer so that our flow reading reflects
  // the current status of the instrument

  MaxCamSerial::Flush(rs232);
  usleep(500000);
  string flowstr=GetFeedback("Flow",5);
  if(debug)  cout << "flowstr=" << flowstr << endl;
  
  return atof(flowstr.c_str());
}

void
MaxCamMFC::setFlowPoint(float ff) {
  // SHOULD ADD SOMETHING THAT KEEPS YOU FROM GOING OVER THE MAXIMUM POSSIBLE
  // SinvValueCCR

  char cmd[32];
  if(ff<10) sprintf(cmd,"Sinv%4.3f",ff);
  if((ff>=10)&&(ff<100)) sprintf(cmd,"Sinv%4.2f",ff);
  if((ff>=100)&&(ff<1000)) sprintf(cmd,"Sinv%4.1f",ff);

  //  append CRC to end of string
  Int_t nbytes=sprintf(cmd,"%s%s",cmd,CalcCRC(cmd).c_str());
  //  cout << "cmd=" << cmd << endl;
  
  MaxCamSerial::Flush(rs232);
  MySerialWrite(rs232,cmd, nbytes);
}


string 
MaxCamMFC::GetFeedback(char* keyword, int arglength)
{
  // read full buffer
  int bytesInInputBuffer=0;
  ioctl(rs232, FIONREAD, &bytesInInputBuffer);
  if(debug)  cout << "bytesInInputBuffer=" << bytesInInputBuffer << endl;
  const int nbytes=bytesInInputBuffer;
  char bufferstr[nbytes];
  MaxCamSerial::SerialRead(rs232, bufferstr, bytesInInputBuffer);
  
  if(debug) printAsReadable(bufferstr,bytesInInputBuffer);

  // sometimes, the command sent to the MFC
  // starts with a '?', but the command sent 
  // back doesn't get the '?' appended to it, 
  // for some reason.  Easy way to deal with 
  // this is to strip '?' off of our search for the response
  // it'll still find the response.
  if(debug) cout << "pre ? strip keyword=" << keyword << endl;
  if(keyword[0]=='?') ++keyword;
  if(debug) cout << "post ? strip keeyord=" << keyword << endl;

  // find the first occurence of the keyword in the char string
  string feedback="";
  char * pch = strstr(bufferstr,keyword);
  // let's assume it's found for now.  We've got to 
  // be able to see at least the keyword.
  if(pch!=NULL){
    
    
      string stringsofar="";
      //    cout << "pch=" << pch << endl;
      unsigned int index=0;
      do{
	if(index>(strlen(keyword)-1)){
	  feedback.push_back(pch[index]);
	  //	cout << "pch[" << index << "]=" << pch[index] << endl;
	  //	cout << "index=" << index << endl;
	  // we're starting to read stuff of interest.
	}

	// my attempt to do my own checksum checking on input
	/*
	// need just the CRC of the handle so far.
	stringsofar.push_back(pch[index]);
	string checksumsofar;
	checksumsofar.push_back(pch[(index+1)]);
	checksumsofar.push_back(pch[(index+2)]);
	checksumsofar.push_back(0x0d);
	checksumsofar.push_back(0x00);
	//      cout << "stringsofar=" << stringsofar << endl;
	//      cout << "checksumsofar=" << checksumsofar.c_str() << endl;
	//      cout << "pch[(index+1)]pch[(index+2)]=" << pch[(index+1)] << pch[(index+2)] << endl;
	//      cout << "CalcCRC(stringsofar.c_str())=" << CalcCRC(stringsofar.c_str()).c_str() << endl;
	if(strcmp(checksumsofar.c_str(),CalcCRC(stringsofar.c_str()).c_str())==0) break;
	*/
	
	++index;
      } while (index<strlen(keyword)+arglength);
      // we've cycled up to the point in the output where 
      // we're reading the response.  syphon it in.
  }
  
  return feedback;
}

string
MaxCamMFC::CalcCRC(const char* c) {
  /*
    This code was taken directly from a Sierra Instruments 
    sample program.  No effort was made to understand it.
    
    Multiplying by 0x0100 just shifts the bits left by two slots
    in other words, 'a' = 0x61
    and so 'a'*0x0100 = 0x6100
    
    ^ is the xor operator.  1^1=0; 1^0=1; 0^0=0.
    3^5 = 6
    0xffff^0x6100 = 0x93ff
    
    & is the and operator.  1&1=1; 1&0=0; 0&0=0;
  */
  
  unsigned int crc;
  crc = 0xffff;
  int ii = 0;
  while (! (c[ii]==0) ) {
    //if(debug) cout << ii << ": " << hex << showbase << crc << dec << noshowbase << endl;
    crc = crc^((unsigned int)(c[ii]*0x0100));
    for (int jj=0; jj<8; jj++) {
      if ( (crc & 0x8000) == 0x8000) {
	crc = ((crc<<1)^0x1021) & 0xffff; 
      } else {
	crc = (crc<<1) & 0xffff;
      }
    }
    ii++;
  }
  //if(debug) cout << "hex:    " << hex << showbase << crc << dec << noshowbase << endl;
  if ( (crc & 0xff00) == 0x0d00) crc += 0x0100;
  if ( (crc & 0x00ff) == 0x000d) crc += 0x0001;
  if ( (crc & 0xff00) == 0x0000) crc += 0x0100;
  if ( (crc & 0x00ff) == 0x0000) crc += 0x0001;
  if(debug) cout << "hex:    " << hex << showbase << crc << dec << noshowbase << endl;
  if(debug) cout << "decimal: " << crc << endl;

  // the old return
  //  return crc;

  // need to output the checksum in a way that's useful.  
  // do it as a string to tack on.
  char checksum0[32];
  // the 04 ensures that later on the checksum
  // gets calculated correctly by padding the number
  // with zeroes at the beginning if necessary!
  sprintf(checksum0,"%04x",crc);
  if(debug) cout << "checksum0=" << checksum0 << endl;
  char checksum1[4];
  checksum1[0]=checksum0[0]; 
  checksum1[1]=checksum0[1];
  char checksum2[4]; 
  checksum2[0]=checksum0[2];
  checksum2[1]=checksum0[3];
  
  if(debug) cout << "checksum1[0]=" << checksum1[0] << endl;
  if(debug) cout << "checksum1[1]=" << checksum1[1] << endl;
  
  if(debug) cout << "checksum2[0]=" << checksum2[0] << endl;
  if(debug) cout << "checksum2[1]=" << checksum2[1] << endl;
  
  // have the hex, need the decimal
  int checksum1decimal, checksum2decimal;
  sscanf(checksum1, "%x", &checksum1decimal); 
  sscanf(checksum2, "%x", &checksum2decimal); 
  
  if(debug) cout << "checksum1decimal=" << checksum1decimal << endl;
  if(debug) cout << "checksum2decimal=" << checksum2decimal << endl;

  string checksum;
  checksum.push_back(checksum1decimal);
  checksum.push_back(checksum2decimal);
  checksum.push_back(0x0d);
  checksum.push_back(0x00);
  //  sprintf(checksum,"%c%c%c%c",checksum1decimal,checksum2decimal,0x0d,0x00);
  if(debug) cout << "checksum=" << checksum << endl;

  return checksum;
}

bool
MaxCamMFC::getDebug(bool b){
  return debug;
}

void
MaxCamMFC::setDebug(bool b){
  debug=b;
}

void 
MaxCamMFC::printBuffer(void){
  // read full buffer
  int bytesInInputBuffer=0;
  ioctl(rs232, FIONREAD, &bytesInInputBuffer);
  cout << "bytesInInputBuffer=" << bytesInInputBuffer << endl;
  const int nbytes=bytesInInputBuffer;
  char bufferstr[nbytes];
  MaxCamSerial::SerialRead(rs232, bufferstr, bytesInInputBuffer);
  
  // the entire buffer is in bufferstr; have to do something about the
  // unprintable characters though.  Loop through and replace them with
  // '@' signs for now.
  for(int i=0; i<bytesInInputBuffer; ++i){
    if((int(bufferstr[i])<33)||(int(bufferstr[i])>126)) bufferstr[i]='@';
  }
  cout << bytesInInputBuffer << endl;
  cout << bufferstr << endl;
}

void 
MaxCamMFC::printAsReadable(char* c, int length){
  for(int i=0; i<length; ++i){
    if((int(c[i])<33)||(int(c[i])>126)) c[i]='@';
  }
  cout << c << endl;
}

void MaxCamMFC::MySerialWrite(SerialHandle rs232, char *sendstring, int nbytes)
{
  MaxCamSerial::SerialWrite(rs232,sendstring,nbytes);
}
