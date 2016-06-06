#include "../../../MaxCam/DmtpcEventTable.hh"
#include "TFile.h"
#include <iostream>

void print_usage()
{
  
  std::cout << "combineRates outputRates [inputRates1] [inputRates2] .. [inputRatesN] " << std::endl;


}

int main(int nargs, char**argv)
{
  
  if (nargs <2) 
  {
    print_usage();
    return 1; 
  }

  TFile * out = new TFile(argv[1], "UPDATE"); 
  DmtpcEventTable * rates, * rates2, * rates3; 
  if (out->Get("tab1") !=NULL)
  {
    rates = (DmtpcEventTable *) out->Get("tab1"); 
  }
  else
  {
    rates = new DmtpcEventTable("tab1","tab1"); 
  }
  if (out->Get("tab2") !=NULL)
  {
     rates2 = (DmtpcEventTable*) out->Get("tab2");  
  }
  else
  {
    rates2 = new DmtpcEventTable("tab2","tab2"); 
  }

  if (out->Get("tab3") !=NULL)
  {
     rates3 = (DmtpcEventTable*) out->Get("tab3");  
  }
  else
  {
    rates3 = new DmtpcEventTable("tab3","tab3"); 
  }
 

  for (int i = 2; i < nargs; i++)
  {
    
    std::cout << " Merging with " << argv[i] << std::endl; 
    TFile f(argv[i]); 
    if (f.Get("tab1")!=NULL)
    {
      rates->combineWith((DmtpcEventTable*) f.Get("tab1")); 
    }
    if (f.Get("tab2")!=NULL)
    {
      rates2->combineWith((DmtpcEventTable*) f.Get("tab2")); 
    }
    if (f.Get("tab3")!=NULL)
    {
      rates3->combineWith((DmtpcEventTable*) f.Get("tab3")); 
    }
  }
   
  out->cd();
  rates->Write();
  rates2->Write();
  rates3->Write();
  out->Write();
  out->Close(); 
}

