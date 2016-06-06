#ifndef CALIB_TOOLS_HH
#define CALIB_TOOLS_HH

#include "../../../MaxCam/DmtpcSkimDataset.hh"
#include "../../../MaxCam/MaxCamClusterImage.hh"
#include <vector>


namespace calibTools
{
   TObjArray* sumImages(DmtpcSkimDataset& ds,TString type); //return an array of the images
   TObjArray* defineRegions(TObjArray* sumImages); //return an array of "masks" for the images
   bool isInRegion(TH2F* region, vector<int> cluster);
   TH2F* makeSpacerMask(TH2F* image, int nspacers, double* m, double* b, double* width);
   TH2F* blurWithoutSpacers(TH2F* image, TH2F* spacermask, double blurfrac);
};


#endif
