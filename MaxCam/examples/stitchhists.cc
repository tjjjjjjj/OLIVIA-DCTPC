#include "TFile.h"
#include "../Dmtpc4ShooterStitcher.hh"
#include "../DmtpcLensCorrection.hh"
#include <vector>
#include <iostream>
#include <TH2.h>
#include "TApplication.h"

int main(int nargs, char **args)
{

  std::vector<const TH2*> images; 

  if (nargs < 3) 
  {
    std::cout << "Usage: ./stitchhists hists.root stitch.root [out.root=hists_stitched.root]"; 
    return 1; 
  }
  

  TFile sfile(args[2]); 
  Dmtpc4ShooterStitcher * stitch = (Dmtpc4ShooterStitcher*) sfile.Get("stitch"); 

  TFile infile(args[1]); 

  for (int i = 0; i < 4; i++)
  {
    images.push_back((const TH2*) infile.Get(stitch->getSerial(i))); 
  }


  TString outname;
   

  if (nargs > 3)
  {
    outname = args[3]; 
  }
  else
  {
    outname = args[1]; 
    outname.ReplaceAll(".root","_stitched.root"); 
  }
  
  TFile out(outname,"RECREATE"); 

  TH2 * result = stitch->stitch(&images); 
  result->Write(); 
  out.Close(); 



}
