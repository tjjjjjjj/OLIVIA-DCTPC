#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/MaxCamClusterImage.hh"
#include <iostream>


int main(int nargs, char ** args)
{

  DmtpcSkimDataset d1; 
  DmtpcSkimDataset d2; 

  d1.openRootFile(args[1]);
  d2.openRootFile(args[2]);


  d1.tree()->SetBranchStatus("*cluster*",0);
  d2.tree()->SetBranchStatus("*cluster*",0);
  d1.tree()->SetBranchStatus("*trigger*",0);
  d2.tree()->SetBranchStatus("*trigger*",0);

  for (int i = 0; i < d1.nevents(); i++)
  {
    d1.getEvent(i);
    d2.getEvent(i);

    for (int c = 0; c < 2; c++)
    {
      if (d1.event()->ntracks(c) != d2.event()->ntracks(c))
      {
        std::cout << "Event " << i << " cam " << c << " differ." << std::endl; 
      }
    }
  }
}
