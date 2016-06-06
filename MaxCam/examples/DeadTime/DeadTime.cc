#include "../../DmtpcDataset.hh"
#include "../../DmtpcEvent.hh"
#include <iostream>
#include <cmath>
#include "../../MaxCamConfig.hh"
#include "TTimeStamp.h"
#include "TChain.h"

using namespace std; 

int main(int nargs, char ** args)
{
  if (nargs < 2)
  {
    cerr << "Need at least one file to calculate" << endl; 
    return 1; 
  }


  for (int i=1; i < nargs; i++)
  {
    double sum = 0; 
    double square_sum = 0; 
    double max = -1;
    double min = 1e20; //big number

    DmtpcDataset ds; 
    ds.openRootFile(args[i]); 
    
    int last_sec = 0; 
    int last_nsec = 0; 
    int n = ds.chain()->GetEntries(); 

    for (int j = 0; j < n; j++)
    {

      ds.getEvent(j); 

      if (j > 1) 
      {
        TTimeStamp* current = ds.event()->UTCtimeStamp() ;

        //cout << current->GetSec() << endl; 

        int secs_diff = current->GetSec() - last_sec;
        int nanosecs_diff = current->GetNanoSec() - last_nsec;
        
        double exposure = 0.001 * ds.event()->ccdConfig(0)->exposureTime; 
        double diff = secs_diff + ((double)nanosecs_diff)/(1e9) - exposure; 

  //      cout << nanosecs_diff<<endl;; 
  //        cout << diff << endl; 
        sum += diff;
        square_sum +=diff*diff; 
        if (diff > max) max = diff;
        if (diff < min) min = diff; 
      }

      last_sec = ds.event()->UTCtimeStamp()->GetSec();
      last_nsec = ds.event()->UTCtimeStamp()->GetNanoSec();

    }
  
    double mean = sum/(n-1); 
    double stdev = sqrt(square_sum/(n-1) - mean*mean); 
    cout << "Stats for \"" << args[i] << "\" (all times in secs)" << endl; 
    cout << "  mean: " << mean << endl; 
    cout << "  stdev: " << stdev << endl; 
    cout << "  max: " << max << endl; 
    cout << "  min: " << min << endl; 
    cout << endl; 

  }

}
