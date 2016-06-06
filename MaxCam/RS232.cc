#include "RS232.hh"
#include <iostream>
#include <fcntl.h>

using std::cout;
using std::endl;
using std::string;

ClassImp(RS232)

RS232::RS232() {}

SerialHandle RS232::Open(string device) 
{
  SerialHandle rs232;
  rs232 = open(device.c_str(), O_RDWR | O_NOCTTY);
  return rs232;
}

int RS232::SetBaudRate(struct termios *rs232_attr, long baudrate) 
{
  switch (baudrate) {
  case 19200:
    rs232_attr->c_cflag |= B19200; break;
  default:
    cout << "unsupported baud rate" << endl;
    return -1;
  }
  return 0;
}

int RS232::Close(SerialHandle rs232) 
{
  close(rs232);
  return true;
}

