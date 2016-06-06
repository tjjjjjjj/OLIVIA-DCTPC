
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
#include "TCanvas.h"

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

//This script needs to open up a list of files, draw out the data trees,
//and reset the scales on the images to match the gain data.


int calibrate(TString rawdatafiles="files.txt",TString keyfilename = "keys.txt",TString calibdir = "./calib/")
{
   TString key = "calib";

   const int nreqmod = 1;
   TString reqmod[nreqmod]={"clean"};

   bool pass;
   DmtpcKeys k(keyfilename,rawdatafiles,key,calibdir,nreqmod,reqmod,pass);
   if(!pass) return -1;

   TFile* routfile=0;
   TTree* calibtree=0;

   for(int f=0; f<k.getNFiles(); f++)
   {
      cout << k.getFile(f) << endl;
      
      //Create DMPTC Database and draw out tree
      DmtpcDataset d;
      d.openRootFile(k.getRootDirName()+k.getFile(f));
      TTree* rawtree = d.tree();

      //Add friends
      k.addFriends(rawtree,f);
      
      //name outputfile, open up outputfile
      TString routfilename = k.getFile(f);
      routfilename.ReplaceAll(".root",key+".root");     
      routfile = new TFile(calibdir+routfilename,"CREATE");
      
      if(routfile->IsOpen())
      {
	 routfile->cd();
      }
      else
      {
	 cout << "Output file already exists; Aborting!" << endl;
	 return -1;
      }

      d.getEvent(0);
      const int ncamera = d.event()->ccdData()->GetEntries();

      //create tree to store calibrated images
      calibtree = new TTree(key,"Calibrated Images");
      TClonesArray* calibimage = new TClonesArray("TH2F",ncamera);

      //create branches
      calibtree->Branch("calibimage","TClonesArray",&calibimage,128000,0);
      gROOT->cd();


      //find old branches
      TClonesArray* cleanimage =0;
      rawtree->SetBranchAddress("cleanimage",&cleanimage);


      //calibration should go here; for now, just reproducing images
      for(int i = 0; i<rawtree->GetEntries(); i++)
//      for(int i = 0; i<50; i++)
      {
	 cout << i << endl;
	 d.getEvent(i);
	 rawtree->GetEntry(i);
	 for(int u=0; u<ncamera; u++)
	 {
	    TH2F* tempimg = (TH2F*)cleanimage->At(u);//->Clone("calibimg");
	    new((*calibimage)[u]) TH2F(*tempimg);
	 }
	 routfile->cd();
	 calibtree->Fill();
	 gROOT->cd();
      }
      routfile->cd();
      calibtree->Write();
      
      delete calibtree;
      routfile->Close();

   }

   return 0;

}
