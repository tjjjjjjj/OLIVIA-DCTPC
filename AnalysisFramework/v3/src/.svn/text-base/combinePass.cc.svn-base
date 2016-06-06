#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "../../../MaxCam/DmtpcSkimEvent.hh"
#include "recoilEnergy.hh"
#include <iostream>
#include "TH2F.h"

double min_Er = 80; 
double max_Er = 200; 

int main(int nargs, char ** args)
{

 
  TFile f(args[1],"RECREATE"); 

  TTree * out = new TTree("skim","skim"); 
  DmtpcSkimEvent * ev = 0;//new DmtpcSkimEvent; 
  DmtpcEvent * old_ev = 0; 
  int run, cam, track; 
  
 // out->SetBranchStatus("*cluster*",0); 
  out->Branch("event","DmtpcSkimEvent",&ev); 
  out->Branch("cam",&cam,"cam/I"); 
  out->Branch("track",&track,"track/I"); 
  out->Branch("run",&run,"run/I"); 
  out->Branch("orig_event","DmtpcEvent",&old_ev); 


  for (int i = 2; i < nargs; i++)
  {
      TFile g(args[i]); 
      
      TTree * in = (TTree*) g.Get("skim"); 
      if (in != 0) 
      {
        std::cout << " File " << args[i] << " contains pass events." << std::endl; 
        ev = 0; 
        old_ev = 0; 
        in->SetBranchAddress("event",&ev);
        in->SetBranchAddress("cam",&cam);
        in->SetBranchAddress("track",&track); 
        in->SetBranchAddress("run",&run); 
        in->SetBranchAddress("orig_event",&old_ev); 
        
        for (int j = 0; j < in->GetEntries(); j++)  
        {
          in->GetEntry(j); 
          double Er = RecoilEnergy::getRecoilEnergy(ev->E(cam,track),cam); 
          std::cout << Er << " keV" << std::endl; 
          if (Er > min_Er && Er < max_Er)
          {
            f.cd(); 
            out->Fill();  
            out->Write(); 
          }
        }
      } 
      g.Close(); 

  }
  f.Close();

  return 0; 
}

