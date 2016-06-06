#include "DmtpcStitchGUIFrame.hh"
#include "../DmtpcRootTools.hh"
#include "TStyle.h"

#include "TApplication.h"


int main (int nargs, char ** args)
{

  TApplication app("app",0,0); 
  DmtpcRootTools::setColorStandard1(); 
//  gDebug = 6; 

  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);

  DmtpcStitchGUIFrame * f = new DmtpcStitchGUIFrame(gClient->GetRoot(),800,800); 


  app.Run();



  return 0; 

}
