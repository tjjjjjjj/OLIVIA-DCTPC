#ifndef RS232_HH
#define RS232_HH

#ifndef __CINT__
#include <termios.h>
#else
struct termios;
#endif

#include <string> 
#include "TNamed.h"

using std::string;

typedef int SerialHandle;

class RS232 : public TNamed {

public:
  // Ctors
  RS232(); 

  static SerialHandle Open(string device); 
  static int SetBaudRate(struct termios *rs232_attr, long baudrate);
  static int Close(SerialHandle rs232);

private:
  ClassDef(RS232,1)
};

#endif

