
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
//and do a clusterfinding and quick reconstruction on them.


int muonRecon(TString rawdatafiles="files.txt", TString keyfilename="keys.txt",TString muondir = "./muon/")
{
   
   TString key = "muon";

   const int nreqmod = 2;
   TString reqmod[nreqmod]={"clean","calib"};

   bool pass;
   DmtpcKeys k(keyfilename,rawdatafiles,key,muondir,nreqmod,reqmod,pass);
   if(!pass) return -1;


   TFile* routfile = 0;
   TTree* muontree = 0;
   for(int f=0; f<keys.getNFiles(); f++)
   {
      //Create DMPTC Database and draw out tree
      DmtpcDataset d;
      d.openRootFile(k.getRootDirName()+k.getFile(f));
      TTree* rawtree = d.tree();

      //Add friends
      k.addFriends(rawtree,f);
      
      //name outputfile, open up outputfile
      TString routfilename = k.getFile(f);
      routfilename.ReplaceAll(".root",key+".root");
      routfile = new TFile(muondir+routfilename,"CREATE");
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
      int ncam = ncamera;

      //create tree to store reconstruction variables
      muontree = new TTree(key,"Reconstruction");
      double range[ncamera][15],x[ncamera][15],y[ncamera][15],E[ncamera][15]; 
      int ntracks[ncamera];


      //Add branches
      muontree->Branch("ncamera",&ncam,"ncamera/I");
      muontree->Branch("E",E,"E[ncamera][15]/D");
      muontree->Branch("range",range,"range[ncamera][15]/D");
      muontree->Branch("x",x,"x[ncamera][15]/D");
      muontree->Branch("y",y,"y[ncamera][15]/D");
      muontree->Branch("ntracks",ntracks,"ntracks[ncamera]/I");
      gROOT->cd();

      TClonesArray* calibimage = 0;
      rawtree->SetBranchAddress("calibimage",&calibimage);

      MaxCamCluster* clust[15];

      TCanvas* c  = new TCanvas("c","c",0,0,1000,500);
      c->Divide(2,1);

      //event selection and reconstruction should go here
//      for(int i = 0; i<rawtree->GetEntries(); i++)
     for(int i = 0; i<50; i++) //for testing
      {
	 cout << i << endl;
	 d.getEvent(i);
	 rawtree->GetEvent(i);
	 for(int u =0; u<ncamera; u++)
	 {
	    //cloned image
	    TH2F* baseimage = (TH2F*)calibimage->At(u)->Clone("baseimage");
	    //rebin, blur, find cluster
	    baseimage->Rebin2D(2,2);
	    baseimage = MaxCamImageTools::blur(baseimage,1,0.8);
	    ntracks[u]=MaxCamImageTools::findClusters(baseimage,
						      clust,120,50000,3,1600);

	    TH2F* clustimg = (TH2F*)baseimage->Clone("clustimg");

	    for(int v =0; v<ntracks[u]; v++)
	    {
	       c->cd(1);
	       baseimage->Draw("colz");
	       vector<int> bins = clust[v]->getCluster();
	       for(int k=0; k<int(bins.size()); k++)
	       {
		  clustimg->SetBinContent(bins[k],100*(v+1));
	       }
	       c->cd(2);
	       clustimg->Draw("colz");		
	       c->Update();

	       //get reconstruction variables
	       double phi = clust[v]->getPhi2();
	       E[u][v] = clust[v]->getIntegral();
	       range[u][v] = clust[v]->getLength2(phi,1);
	       clust[v]->getXY(x[u][v],y[u][v]);

	    }
	 }
	 routfile->cd();
	 muontree->Fill();
	 gROOT->cd();
      }
     routfile->cd();
     muontree->Write();
     
     delete muontree;
     routfile->Close();
     
   }

   return 0;
   
}
