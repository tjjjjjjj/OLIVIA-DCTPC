#include "../../DmtpcSkimDataset.hh"
#include <iostream>

int main (int nargs, char ** args) 
{

  int last_run_number = -1; 
  DmtpcSkimDataset d; 
  
  for (int i = 1; i < nargs; i++)
  {
    std::cout << "Loading dataset " << args[i] << std::endl;
    d.openRootFile(args[i]); 
    d.getEvent(0); 
    if (last_run_number == d.event()->runNumber())
    {
      std::cerr << "REASON: RUN NUMBER DIDN'T INCREMENT " << std::endl;  
      return 1; 
    }
  
    last_run_number = d.event()->runNumber(); 
  }

  return 0; 

}
