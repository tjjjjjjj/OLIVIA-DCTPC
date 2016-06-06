#include "cleanSkimFunctions.hh"
#include <cstdlib>
#include <vector>
#include "TApplication.h" 
#include "TCanvas.h"
#include "DmtpcRootTools.hh" 

int main(int nargs, char ** args)
{

  TApplication app("app",0,0); 
  int run;
  if (nargs > 2) run = atoi(args[2]); 


  CleanSkimConfig cfg(args[1]);

  DmtpcGainMap * gm; 
  Dmtpc4ShooterStitcher * stitch = cleanSkimFunctions::loadStitch(run,&cfg); 
  vector<string> ids; 
  ids.push_back("A80333"); 
  ids.push_back("A80334"); 
  ids.push_back("110121"); 
  ids.push_back("100534"); 
  cleanSkimFunctions::loadGainMaps(&gm,0,4,&cfg, stitch, &(ids[0])); 

  DmtpcRootTools::setColorStandard1(); 
  TCanvas c; 
  c.cd(); 
  gm->getGainMap()->Draw("colz"); 
  app.Run(); 


  return 0; 

}

