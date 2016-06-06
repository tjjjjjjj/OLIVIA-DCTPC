#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include <iostream>
#include <vector>

int main (int nargs, char ** args)
{

  if (nargs < 2)  
  {
    std::cout << "Give me at least one file please. " << std::endl; 
  }

  for (int arg = 1; arg < nargs; arg++)
  {

    DmtpcSkimDataset d; 
    d.openRootFile(args[arg]); 
    d.tree()->SetBranchStatus("*cluster*",0); 
    d.tree()->SetBranchStatus("*trigger*",0); 

    d.getEvent(0); 

    std::vector<int> cam_counts(d.event()->ncamera(),0); 
    
    std::vector<std::vector<int> > cam_distributions(d.event()->ncamera(), std::vector<int>(15,0)) ; 

    for (int i = 0; i < d.nevents(); i++)
    {
      d.getEvent(i);
      for (int c = 0; c < d.event()->ncamera(); c++)
      {
        cam_counts[c] += d.event()->ntracks(c); 
        cam_distributions[c][d.event()->ntracks(c)]++;
      }
    }

    std::cout << "File: " << args[arg] << std::endl; 
    for (int c = 0; c < d.event()->ncamera(); c++)
    {
      std::cout << "\tCamera " << c << ": " << cam_counts[c] << " Total Tracks. " << std::endl; 
      std::cout << "\t\tDistribution:" << std::endl; 
      for (int t = 0; t < 15; t++)
      {
        if (cam_distributions[c][t]>0)
        {
          std::cout << "\t\t\t" << t << " Clusters: " << cam_distributions[c][t] <<std::endl;
        }
      }
    }

  }
}
