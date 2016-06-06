#ifndef MAXCAM_MFC_HH
#define MAXCAM_MFC_HH

#include <string> 
#include <iostream>
#include "TNamed.h"
#include "MaxCamSerial.hh"

using std::cout;
using std::endl;
using std::string;
using std::hex;
using std::dec;
using std::showbase;
using std::noshowbase;

// A class for a Mass Flow Controller 
// e.g. the Sierra Instruments C100

class MaxCamMFC : public TNamed {

public:

  // Ctors
  MaxCamMFC() {};

  MaxCamMFC(const char *name, const char *title, char *port="/dev/ttyS0", string model="SIERRAC100", char *dbaccess=0);

  MaxCamMFC(const MaxCamMFC &other);

  virtual ~MaxCamMFC() { cout << "closing serial port" << endl; MaxCamSerial::Close(rs232); };

  double flowRate() { return currentFlowRate; }
  double currentFlowRate;

  long getBaudRate() { return baudRate; }

  void openValve(void);
  void closeValve(void);
  void purgeValve(void);
  void setValve(int ii);
  void setUnits(int ii);
  void setGas(int ii);
  void setSource(int ii);
  void setOutput(int ii);
  void setFlowPoint(float ff);

  // set internal variables
  void setPause(float p);
  
  float getFlow(void);
  string queryMFC(char* c,int arglength);
  void printBuffer(void);
  void printAsReadable(char* c,int length);

  string CalcCRC(const char* c);
  string getVersion(void);

  int getSerial(void);
  float getFactor(void);
  string getFlowPoint(void);
  int getSource(void);
  string getPassword(void);
  int getValve(void);
  int getOutput(void);
  int getUnits(void);
  float getMax(void);
  int getGas(void);
  
  string GetFeedback(char* keyword, int arglength);
  
  // get commands
  bool getDebug(bool b);

  // set commands
  void setDebug(bool b);

  void Fill(float rate, int seconds);

  void MySerialWrite(SerialHandle rs232, char *sendstring, int nbytes);
  
  
private:

  //MaxCamSerial mcs;

  int baudRate;
  string port;
  bool debug;
  float pause; // how much to pause on send commands

  SerialHandle rs232;

  ClassDef(MaxCamMFC,1)
};

#endif

