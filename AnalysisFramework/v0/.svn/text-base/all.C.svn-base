
#include "../../MaxCam/MaxCamMC.hh"
#include "../../MaxCam/MaxCamImageTools.hh"
#include "../../MaxCam/MaxCamCluster.hh"
#include "../../MaxCam/DmtpcEvent.hh"
#include "../../MaxCam/DmtpcDataset.hh"
#include "../../MaxCam/DmtpcKeys.hh"
#include "TF1.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"
#include "TRandom3.h"

#include <iostream>
#include <vector>

using namespace std;

int all(TString rawdatafiles = "files.txt",TString keyfilename ="keys.txt",
	TString cleandir = "./clean/",TString calibdir = "./calib/", 
	TString recondir = "./recon/")
{

   
   // Load up the scripts
   gROOT->ProcessLine(".L cleanData.C+");
   gROOT->ProcessLine(".L calibrate.C+");
   gROOT->ProcessLine(".L quickRecon.C+");
//   gROOT->ProcessLine(".L cuts.C");

   int clean = cleanData(rawdatafiles,keyfilename,cleandir);
   TString cont;
   if (clean!=0)
   {
      cout << "cleanData failed! Attempt to continue, Y/N?";
      cin >> cont;
   }

   if(clean == 0 || cont == "Y" || cont == "y")
   {
      int calib = calibrate(rawdatafiles,keyfilename,calibdir);
   }
   else
   {
      cout << "Exiting!" << endl;
      return -1;
   }

   if(calib != 0)
   {
      cont = "";
      cout << "calibrate failed! Attempt to continue, Y/N?";
      cin >> cont;
   }

   if(calib == 0 || cont =="Y" || cont == "y")
   {
      int recon = quickRecon(rawdatafiles,keyfilename,recondir);
   }
   else
   {
      cout << "Exiting!" << endl;
      return -1;
   }

   if(recon != 0)
   {
      cont = "";
      cout << "quickRecon failed! Attempt to continue, Y/N?";
      cin >> cont;
   }


   if(recon == 0 && calib == 0 && clean == 0)
   {
      cout << "Sucessful completion of loaded scripts!" << endl;
      return 0;
   }
   else
   {
      cout << "Exiting with errors!" << endl;
      return -1;
   }
   

}
