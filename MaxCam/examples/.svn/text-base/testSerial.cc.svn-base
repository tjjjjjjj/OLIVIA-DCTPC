// test serial communications with the Sierra Instruments Mass Flow Controller
// (model = C100)

#include <iostream>
#include <fcntl.h>
#include <termios.h>
#include <unistd.h>

#include <ctype.h>  // isdigit()
#include <string>

using std::cout;
using std::endl;
using std::string;

// set-up for the Sierra Instruments MFC
#define BAUDRATE B9600
#define DEVICE "/dev/ttyS0"

int main() {

  int fd;

  //----------------------------------------------------------
  fd = open(DEVICE, O_RDWR | O_NOCTTY | O_NDELAY);
  //fd = open(DEVICE, O_RDWR | O_NOCTTY);
  if (fd < 0) {
    cout << ": Could not connect to " << DEVICE << endl;
    return fd;
  }
  else 
    cout << ": Connected to " << DEVICE << endl;
  //----------------------------------------------------------

  struct termios options;
  if (tcgetattr(fd, &options)) {
    cout << ": tcgetattr failed " << endl;
    return -1;
  }// else {
   // cout << "tcgetattr succeeded" << endl;
  // }

  /* 
     BAUDRATE: Set bps rate. You could also use cfsetispeed and cfsetospeed.
     CRTSCTS : hardware handshaking (only used if the cable has all 
     necessary lines -- the Sierra C100 does not)
     CS8     : 8n1 (8bit,no parity,1 stopbit)
     CLOCAL  : local connection, no modem contol
     CREAD   : enable receiving characters
  */
  options.c_cflag = BAUDRATE | CS8 | CLOCAL | CREAD; 

  // behaves like a 1 second timeout
  options.c_cc[VMIN]  = 0;   
  options.c_cc[VTIME] = 10;  // t=TIME*0.1s, so 10 is 1 second

  cfsetospeed(&options, BAUDRATE);
  cfsetispeed(&options, BAUDRATE);

  // Clear the line
  tcflush(fd,TCIFLUSH);
  
  // Set the new options for the port...
  if (tcsetattr(fd, TCSANOW, &options)) {
    cout << "tcsetattr failed\n" << endl;
    close(fd);
    return -1;
  }

  //
  // Terminal settings done, now handle input
  //

  // MFC continually broadcasts Flow values
  // e.g. 'Flow3.111ccr'
  // where cc is a checksum and r is a carriage return

  char c[1024];
  int bytes_read   = 0;
  int nBytesToRead = 15;
  string cmdStr, valStr;
  int index;
  bool validStart;

  while (1) {
    // read() returns the number of bytes *actually* read
    bytes_read = read(fd, c, nBytesToRead);
    // responses begin with a letter or a ?
    validStart = isalpha(c[0]) || (c[0]=='?');
    // if you have a valid response, parse it
    if ( bytes_read > -1 && validStart) {
      cout << bytes_read << ":  [";
      // debugging, just print out the entire response
      for (int ii=0; ii<bytes_read; ii++) {
	cout << c[ii]; 
      }
      cout << "]" << endl;

      index = 0;
      // look for command names, could be letters or question mark
      cmdStr.clear();
      for (int ii=0; ii<bytes_read; ii++) {
  	if (!(isalpha(c[ii]) || (c[ii]=='?')) )
  	  break;
  	index++;
  	cmdStr += c[ii];
      }
      cout << "[" << cmdStr << "]" << endl;

      // Next comes the value 
      // which is a number (can have a decimal point or be negative)
      valStr.clear();
      for (int ii=index; ii<bytes_read; ii++) {
  	if ( ! (isdigit(c[ii]) || c[ii]=='.' || c[ii]=='-') )
  	  break;
  	index++;
  	valStr += c[ii];
      }
      cout << "[" << valStr << "]" << endl;

      //// then read two bytes into checksum
      //for (int ii=0; ii<2; ii++) {
      //	index++;
      //}
      //// then read carriage return
      //char carRet;
      //carRet = c[index];
      //cout << (c[index]=='\n') << (c[index]=='\r') << endl;
      //cout << "[" << carRet << "]" << endl;

      cout << endl;
    }
  }

  //std::string buffer;
  //buffer.clear();
  //char ENQ=5;
  //bytes_send = write(fd, &ENQ, 1 );
  //assert (bytes_send==1);
  //
  //int ncount=0;
  //while (1) {
  //  int nread = read(fd, &c[ncount], 1 );
  //  if (nread<1) continue;
  //  buffer += c[ncount];
  //  if (++ncount==15) break;
  //}
  //
  //int flg;
  //sscanf( buffer.c_str(), "%d, %le", &flg, &currentValue);
  //switch(flg) {
  //  //case 0: cout << GetName() << ": Current value is "     << currentValue << endl; break;
  //case 1: cout << GetName() << ": Underrange, P set to " << currentValue << endl; break;
  //case 2: cout << GetName() << ": Overrange, P set to "  << currentValue << endl; break;
  //}
  //currentValueRMS=0;

  
  close(fd);

  return 0;
}

