#include "../../DmtpcSkimDataset.hh"
#include <iostream>

int main (int nargs, char ** args) 
{

  for (int t = 0; t < atoi(args[2]); t++)
  {
    std::cout << "iteration " << t << std::endl; 
    DmtpcSkimDataset d; 
    d.openRootFile(args[1]); 
    for (int i = 0; i < 100; i++)
    {
      if (i >= d.nevents()) break; 
      d.getEvent(i); 
      assert(d.event("skim")); 
      assert(d.event("ci")); 
    }
    gROOT->Print(); 
  }

  return 0; 

}
