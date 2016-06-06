#include "../MaxCamImageTools.hh"
#include "../DmtpcIterators.hh"
#include "../DmtpcDataset.hh"
#include <vector>
#include <stdio.h>
#include <iostream>
#include "TH2.h"
#include "TFile.h" 


int main(int nargs, char ** args) 
{


  DmtpcDataset  *d = new DmtpcDataset; 
  d->openRootFile(args[1]); 
  int n = 100;
  TFile out (args[2], "RECREATE"); 
  for (int cam = 0; cam < 3; cam++)
  {
    DmtpcDatasetImageIterator beg(d,cam,0); 
    DmtpcDatasetImageIterator end(d,cam,n); 
    TH2 * median = MaxCamImageTools::histStackNthElement<unsigned short> (beg, end, 0,n,65536,true); 
    char buf[20];
    sprintf(buf,"med%d",cam); 
    median->SetName(buf); 
    median->SetTitle(buf); 
    out.cd(); 
    median->Write(); 
    
    sprintf(buf,"bias%d",cam); 
    TH2 * bias = NEW_HIST2D_WITH_SAME_SIZE(median,TH2F,buf); 
    TH2 * bias_rms = NEW_HIST2D_WITH_SAME_SIZE(median,TH2F,buf); 
    for (int i = 0; i < n; i++)
    {
      d->getEvent(i);  
      bias->Add(d->event()->ccdData(cam), 1./n); 
    }
    bias->Write();
    sprintf(buf,"diff%d",cam); 
    TH2 * diff = (TH2*) bias->Clone(buf); 
    diff->SetTitle(buf); 
    diff->Add(median,-1); 
    diff->Write(); 
    out.Flush(); 
    cout << "WROTE" << endl; 
  }

  

  out.Close(); 

}
