#include "TFile.h"
#include "../../../MaxCam/DmtpcEventTable.hh"
#include <iostream>

int main (int nargs, char ** args)
{

  for (int i = 1; i < nargs; i++)
  {
    TFile f(args[i]);
    if (f.Get("tab1")!=NULL)
    {

      std::cout << "Printing table1 from file " << args[i] <<std::endl; 
      ((DmtpcEventTable*)f.Get("tab1"))->print(); 
      std::cout << std::endl << "Printing table2 from file " << args[i] <<std::endl; 
      ((DmtpcEventTable*)f.Get("tab2"))->print(); 
      std::cout << std::endl << "Printing table3 from file " << args[i] <<std::endl; 
      ((DmtpcEventTable*)f.Get("tab3"))->print(); 
    }
  }

  return 0; 

}
