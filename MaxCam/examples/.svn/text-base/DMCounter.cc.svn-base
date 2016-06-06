#include "../LabJackU3.hh"
#include <unistd.h>
#include <iostream>
int main(int nargs, char ** args)
{
  LabJackU3 u3; 

  u3.setUpCounters(8); 	
  double value; 
  while(true) 
  {
    sleep(1); 
    u3.readCounters(&value,0); 
    std::cout << "Counter Value: " << value << std::endl;  

    if (value > 10) 
    {
      u3.resetCounters(true,false);  	

    }
  }

  return 0; 
}	
