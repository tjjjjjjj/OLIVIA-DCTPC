#include "../../DmtpcDataset.hh"
#include "../../DmtpcEvent.hh"
#include <iostream>
#include <cmath>
#include "../../MaxCamConfig.hh"
#include "TDatime.h"
#include "TChain.h"

using namespace std; 

static const int m = 100; 

int main(int nargs, char ** args)
{

  if (nargs < 2)
  {
    cerr << "Need at least one file to calculate" << endl; 
    return 1; 
  }


  for (int i = 1; i < nargs; i++)
  {
    double sum = 0; 
    double c= 0; 

    DmtpcDataset ds; 
    ds.openRootFile(args[i]); 

    time_t last_time = 0; 

    int n = ds.chain()->GetEntries(); 

    for (int j = 0; j < n; j+=m)
    {
      ds.getEvent(j);
      TDatime* current = ds.event()->timeStamp();
      if (j>1)
      {
      //  cout << current->Convert() << endl; 
        double diff = current->Convert() - last_time;   
        diff/=m;
        double exposure = 0.001 * ds.event()->ccdConfig(0)->exposureTime; 
        diff-=exposure;
        sum+=diff; 
        c++; 
      }
      last_time = current->Convert(); 
    }

    double mean = sum/(c-1); 
    cout << "Stats for \"" << args[i] << "\" (in secs, using interval "<<m<<")" << endl;
    cout << "  mean: " << mean << endl; 
    cout <<endl; 
  }

}
