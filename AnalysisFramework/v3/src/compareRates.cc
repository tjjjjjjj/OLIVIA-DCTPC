#include <TString.h>
#include "TFile.h"
#include "../../../MaxCam/DmtpcEventTable.hh"
#include <iostream>

int main(int nargs, char ** args)
{
   TString prefix1(args[1]);
   TString prefix2(args[2]);

  
   for (int i = 801; i < 1032; i++)
   {
      TString fn1 = prefix1;
      fn1+=i;
      fn1+= ".root" ;
      TString fn2 = prefix2;
      fn2+=i;
      fn2+= ".root" ;

      TFile f1(fn1.Data());
      TFile f2(fn2.Data());
  
      DmtpcEventTable * t1 = (DmtpcEventTable *) f1.Get("tab1");
      DmtpcEventTable * t2 = (DmtpcEventTable *) f2.Get("tab1");
      
      
      std::cout <<  (t1->value("poscut","NTomPass_0") - t2->value("poscut","NCosminPass0") + 
          t1->value("poscut","NTomPass_1") - t2->value("poscut","NCosminPass1")) << std::endl; 
   }

   return 0; 
}
