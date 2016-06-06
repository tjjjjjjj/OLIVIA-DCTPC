
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
#include "TROOT.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

//This script needs to open up a list of files, draw out the data histograms,
//evaluate whether any are sparks, subtract the pedestal, and remove lone high pixels.
//The cleaned image will be added to the tree and saved elsewhere


int cleanData(TString rawdatafiles="files.txt",TString keyfilename = "keys.txt", TString cleandir="./clean/")
{

   TString key = "clean";

   const int nreqmod = 0;
   TString reqmod[nreqmod];
   
   bool pass;
   DmtpcKeys k(keyfilename,rawdatafiles,key,cleandir,nreqmod,reqmod,pass);
   if(!pass) return -1;

   TCanvas* c = new TCanvas("c","c",0,0,1000,1000);
   c->Divide(2,2);
   c->cd(1);

   TFile* routfile;
   TTree* cleantree;
   TClonesArray* cleanimage;
   TH2F* tempimg;

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
      const int nbinsx = d.event()->ccdData(0)->GetNbinsX();
      const int nbinsy = d.event()->ccdData(0)->GetNbinsY();

      //create tree to store clean images and sparking data
      cleantree = new TTree(key,"Cleaned Images");
      cleanimage = new TClonesArray("TH2F",ncamera);
      Bool_t spark[ncamera];
      double integral[ncamera];

      //Add branches
      cleantree->Branch("cleanimage","TClonesArray",&cleanimage,128000,0);
      cleantree->Branch("spark",&spark,"spark[2]/O");
      cleantree->Branch("integral",&integral,"integral[2]/D");
      gROOT->cd();

      //open up file for bias frame saving
      TString biasoutfilename = routfilename;
      biasoutfilename.ReplaceAll(key+".root","bias.root");

      TFile* biasoutfile = new TFile(cleandir+biasoutfilename,"RECREATE");

      TH2F* biasframe = d.event()->ccdData(0);
      TTree* biastree = new TTree("bias","Bias Information");

      biastree->Branch("biasframe","TH2F",&biasframe,128000,0);

      gROOT->cd();

      double threshhold[ncamera];
      int nframes[ncamera];
      for(int i=0; i<ncamera; i++){threshhold[i]=200; nframes[i]=0;}

      TH2F* secondarybias[ncamera];

      //find the spark cut
      const int nev = rawtree->GetEntries();
      double countthresh[nev];
      double sparkcut[ncamera];


      for(int j=0; j<ncamera; j++)
      {
	 cout << j << ":" << endl;
	 countthresh[j] =d.getBiasFrame(j+1)->Integral()/(65536);
	 sparkcut[j] =  MaxCamImageTools::countPixelsAboveThreshold(d.getBiasFrame(j+1),countthresh[j])-2000;

	 cout << "count: " << countthresh[j] << "\t cut: " << sparkcut[j] << endl;

	 
      }
      
      //Grab every 100th image to make bias frame; must check for sparks
      for(int i=0; i<nev; i=i+100)
      {
	 cout << i << endl;
	 d.getEvent(i);
	 for(int u=0; u<ncamera; u++)
	 {
	    tempimg = (TH2F*)d.event()->ccdData(u)->Clone("backimg");
//	    MaxCamImageTools::killLonePixels2(tempimg,1300);
	    int nabove = MaxCamImageTools::countPixelsAboveThreshold(tempimg,countthresh[u]);
	    if(nabove<sparkcut[u])
	    {
	       if(nframes[u]==0) 
	       {
		  secondarybias[u]=(TH2F*)tempimg->Clone("average");
	       }
	       else
	       {
		  secondarybias[u]->Add(tempimg,1);
	       }
	       nframes[u]++;
	    }
	    
	 }
	 
      }

      gROOT->Delete("backimg;*");

      for(int i=0; i<ncamera; i++)
      {
	 secondarybias[i]->Scale(1/double(nframes[i]));
 	 biasframe = secondarybias[i];
 	 biasoutfile->cd();
 	 biastree->Fill();
 	 gROOT->cd();
	 
      }

      cout << "All preclean activities done" << endl;
      c->Update();

      //cleaning should go here
      for(int i = 0; i<nev; i++)
//      for(int i = 0; i<100; i++) //For testing
      {
	 cout << i << endl;
	 d.getEvent(i);
	 for(int u=0; u<ncamera; u++)
	 {
//	    c->cd(u+1);

	    
	    //create cloned image
	    tempimg = (TH2F*)d.event()->ccdData(u);

	    //check for sparking
	    int nabove = MaxCamImageTools::countPixelsAboveThreshold(tempimg,countthresh[u]);

	    if(nabove < sparkcut[u])
	       spark[u]=0;
	    else
	       spark[u]=1;

	    if(spark[u] == 1) cout << "SPARK!" << endl;


	    //subtract of temp bias
	    tempimg->Add(secondarybias[u],-1);
	    //kill lone pixels
	    MaxCamImageTools::killLonePixels2(tempimg,threshhold[u]);
   
	    
	    //additionally correct
	    double perpx = tempimg->Integral()/(nbinsx*nbinsy);

	    if(spark[u]==0)
	       MaxCamImageTools::subtractPedestal(tempimg,perpx);

	    integral[u] = tempimg->Integral();
	    
	    //add to TClonesArray
	    new((*cleanimage)[u]) TH2F(*tempimg);
	    


	 }
	 routfile->cd();
	 cleantree->Fill();
	 gROOT->cd();

      }

      routfile->cd();
      cleantree->Write();

      delete cleantree;

      routfile->Close();

      biasoutfile->cd();
      biastree->Write();
      
      delete biastree;
      
      biasoutfile->Close();

   }

      
   return 0;
}
