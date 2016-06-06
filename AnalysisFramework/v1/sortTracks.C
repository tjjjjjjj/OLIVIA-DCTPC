
#include "../../MaxCam/MaxCamMC.hh"
#include "../../MaxCam/MaxCamImageTools.hh"
#include "../../MaxCam/MaxCamCluster.hh"
#include "../../MaxCam/DmtpcEvent.hh"
#include "../../MaxCam/DmtpcDataset.hh"
#include "../../MaxCam/DmtpcKeys.hh"
#include "TF1.h"
#include "TH3.h"
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
#include "TROOT.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

//This script needs to open up a list of files, draw out the data histograms,
//evaluate whether any are sparks, subtract the pedestal, and remove lone high pixels.
//The cleaned image will be added to the tree and saved elsewhere

double isAlpha(MaxCamCluster* track);
double isNeutron(MaxCamCluster* track);
double isWorm(MaxCamCluster* track);


int sortTracks(TString rawdatafiles="files.txt",TString keyfilename = "keys.txt", TString cleandir="./sort/")
{
   //This module runs on skim trees, must have sk in front of key
   TString key = "sksort";

   const int nreqmod = 1;
   TString reqmod[nreqmod] = {"skim"};
   
   bool pass;
   DmtpcKeys k(keyfilename,rawdatafiles,key,cleandir,nreqmod,reqmod,pass);
   if(!pass) return -1;

   TFile* routfile;
   TTree* sorttree;

   for(int f=0; f<k.getNFiles(); f++)
   {
      //Create DMPTC Database and draw out tree
      DmtpcDataset d;
      d.openRootFile(k.getRootDirName()+k.getFile(f));
      TTree* rawtree = k.getBaseTree(f);
      k.addFriends(rawtree,f);

      //name outputfile, open up outputfile
      TString routfilename = k.getFile(f);
      routfilename.ReplaceAll(".root",key+".root");
      routfile = new TFile(cleandir+routfilename,"CREATE");
      
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

      //Get tracks out of the tree
      TObjArray* clustarr=0;
      TClonesArray* cleanimg=0;
      rawtree->SetBranchAddress("clusters",&clustarr);
      rawtree->SetBranchAddress("cleanimage",&cleanimg);

      //create tree to store sorting data
      sorttree = new TTree(key,"Sorting");
      //Each track in each camera needs a place
      double alpha[2][15], frecoil[2][15],worm[2][15];
      double other[2][15];


      //Add branches
      sorttree->Branch("alpha",&alpha,"alpha[2][15]/D");
      sorttree->Branch("frecoil",&frecoil,"frecoil[2][15]/D");
      sorttree->Branch("worm",&worm,"worm[2][15]/D");
      sorttree->Branch("other",&other,"other[2][15]/D");
      gROOT->cd();

//       TCanvas* c0 = new TCanvas("c0","c0",0,0,800,800);
//       c0->Divide(2,2);
//       c0->cd();
      
      TString temp;
      //sorting should go here
      for(int i = 0; i<rawtree->GetEntries(); i++)
//      for(int i = 0; i<100; i++) //For testing
      {
	 cout << i << endl;
	 d.getEvent(i);
	 rawtree->GetEvent(i);
	 TObjArray* c;
	 //loop over cameras
	 for(int j=0; j<clustarr->GetEntries(); j++)
	 {
	    c = (TObjArray*)clustarr->At(j);
	    TH2F* img = (TH2F*)cleanimg->At(j);
	    for(int k=0; k<c->GetEntries(); k++)
	    {
	       MaxCamCluster* track = (MaxCamCluster*)c->At(k);

	       vector <int> px = track->getCluster();
	       cout << "Npx: " << px.size() << endl;

	       alpha[j][k] = isAlpha(track);
	       frecoil[j][k] = isNeutron(track);
	       worm[j][k] = isWorm(track);
	       if(worm[j][k] >0)
		  cout << "WORM!" << endl;

	       if(alpha[j][k]+frecoil[j][k]+worm[j][k] == 0)
	       {cout << "other \n"; other[j][k] = 1;}

	       vector <float> binvals;

	       TH2F* tempimg = (TH2F*)img->Clone();
	       for(int l=0; l<int(px.size()); l++)
	       {
		  binvals.push_back(img->GetBinContent(px[l]));
	       }
	       double maxval =  *max_element(binvals.begin(),binvals.end());
	       cout << maxval << endl;


	    }
	 }
	 routfile->cd();
	 sorttree->Fill();
	 gROOT->cd();

      }



      routfile->cd();
      sorttree->Write();
      delete sorttree;

      routfile->Close();
   }

      
   return 0;
}

double isWorm(MaxCamCluster* track)
{
   cout << "Energy density: " << track->getEnergyDensity() << endl;

   if(track->getEnergyDensity() > 18 || track->getCluster().size() <= 64)
      return 1;
   else
      return 0;
}
double isNeutron(MaxCamCluster* track)
{
   return 0;
}
double isAlpha(MaxCamCluster* track)
{
   cout << "Integral: " << track->getIntegral() << endl;
   
   if(track->getIntegral() > 3000 && track->getEnergyDensity() < 20)
      return 1;
   else
      return 0;
}
