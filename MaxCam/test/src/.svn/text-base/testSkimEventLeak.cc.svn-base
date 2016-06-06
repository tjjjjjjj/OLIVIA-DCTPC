#include "../../DmtpcSkimEvent.hh"
#include <iostream>
#include "TChain.h"
#include "TFile.h"
#include <assert.h>

int main (int nargs, char ** args) 
{

  for (int t = 1; t < atoi(args[2]); t++)
  {
    TFile * f = TFile::Open(args[1]); 
    DmtpcSkimEvent * event = new DmtpcSkimEvent(); 
    DmtpcSkimEvent * event2 = new DmtpcSkimEvent(); 
    TChain * ch = new TChain("skim"); 
    ch->SetBranchAddress("event",&event); 
    ch->Add(args[1]); 
    TChain * ch2 = new TChain("ci"); 
    ch2->SetBranchAddress("event2",&event); 
    ch2->Add(args[1]); 
    for (int i = 0; i < 20; i++)
    {
      if (i >= ch->GetEntries()) break; 
      ch->GetEvent(i); 
      ch2->GetEvent(i); 
      assert(event); 
      assert(event2); 
    }

    delete event; 
    delete event2; 
    delete ch; 
    delete ch2; 
    f->Close();
    delete f; 
  }

  return 0; 

}
