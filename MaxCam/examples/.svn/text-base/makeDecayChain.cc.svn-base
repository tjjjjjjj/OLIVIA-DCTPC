#include "TFile.h"
#include "../DmtpcDecayChain.hh"
#include <stdio.h>


int main (int nargs, char ** args)
{

  if (nargs < 3 ) 
  {
    printf("usage: makeDecayChain in.chain out.root [maxtime=86400] [binning = 10] [dt = 0.01] "); 
    return 0;
  }

  double maxtime = nargs > 3 ? atof(args[3]) : 86400; 
  double binning = nargs > 4 ? atof(args[4]) : 10; 
  double dt = nargs > 5 ? atof(args[5]) : 0.01; 
  TFile out (args[2],"RECREATE"); 

  DmtpcDecayChain *c = new DmtpcDecayChain(args[1],maxtime,binning,dt); 
  out.cd(); 
  c->Write("chain"); 
  out.Close(); 

  return 0; 
}
