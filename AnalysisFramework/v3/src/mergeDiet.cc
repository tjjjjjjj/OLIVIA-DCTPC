#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include <iostream>
#include "TFile.h"
#include <strings.h>
const static char * ext = ".root";
//                      0   1   2   3   4    5          7         9         
int run_boundaries[] = {801,851,901,951,1001,1032,1074, 779, 785, 100421001, 
//                     10                    12                    14 
                       100421011, 100421021, 100422001, 100422005, 100429002, 100429011,
//                     16                    18         19                    21         22 
                       100502001, 100502009, 100429002, 100429006, 100429011, 100502001, 100502005, 100502009
};




int main (int nargs, char ** args)
{
  int n = atoi(args[1]);  
  
  TString name = "merged";
  name+=run_boundaries[n]; 
  name+="_";
  name+=run_boundaries[n+1]-1; 
  name+=".root";

  TFile * f = new TFile(name, "RECREATE"); 
  TTree * skim = new TTree("skim","skim"); 
  DmtpcSkimEvent * ev = 0;  
  skim->Branch("event","DmtpcSkimEvent",&ev); 

  

  for (int i = run_boundaries[n]; i < run_boundaries[n+1];  i++)
  {
    DmtpcSkimDataset ds; 
    char * skimfile = (char *) malloc(strlen(args[2])+40); 
    if (i > 99999) //For MC stuff 
      sprintf(skimfile,"%s/dmtpc_run%dskim.root",args[2],i); 
    else
      sprintf(skimfile,"%s/dmtpc_run%05dskim.root",args[2],i); 
    ds.openRootFile(skimfile);  
    ds.loadClusters(false); 

    std::cout << "Now processing " << skimfile <<std::endl;
    for (int j = 0; j < ds.nevents(); j++)
    {
      if (j % 100 == 1) std::cout << j << std::endl; 
      ds.getEvent(j); 
      ev = new DmtpcSkimEvent(*ds.event(),true); 
      f->cd(); 
      skim->Fill();
    }

    skim->Write(); 
  }

  f->Close();
}

  
