
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


int quickRecon(TString rawdatafiles="files.txt", TString keyfilename="keys.txt",TString recondir = "./recon/")
{
   
   TString key = "recon";

   const int nreqmod = 1;
   TString reqmod[nreqmod]={"clean"};

   bool pass;
   DmtpcKeys k(keyfilename,rawdatafiles,key,recondir,nreqmod,reqmod,pass);
   if(!pass) return -1;

   TFile* routfile;
   TTree* recontree;

   for(int f=0; f<k.getNFiles(); f++)
   {
      //Create DMPTC Database and draw out tree
      DmtpcDataset d;
      d.openRootFile(k.getRootDirName()+k.getFile(f));
      TTree* rawtree = d.tree();

      k.addFriends(rawtree,f);
      
      
      //name outputfile, open up outputfile
      TString routfilename = k.getFile(f);
      routfilename.ReplaceAll(".root",key+".root");
      routfile = new TFile(recondir+routfilename,"CREATE");

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
      recontree = new TTree(key,"Reconstruction");
      double theta[ncamera][15], phi[ncamera][15], E[ncamera][15];
      double range[ncamera][15], x[ncamera][15], y[ncamera][15];
      double skewness[ncamera][15];
      bool edge[ncamera][15];
      int ntracks[ncamera];
      int date;
      int time;


      //Add branches
      recontree->Branch("ncamera",&ncam,"ncamera/I");
      recontree->Branch("theta",theta,"theta[ncamera][15]/D");
      recontree->Branch("phi",phi,"phi[ncamera][15]/D");
      recontree->Branch("E",E,"E[ncamera][15]/D");
      recontree->Branch("range",range,"range[ncamera][15]/D");
      recontree->Branch("x",x,"x[ncamera][15]/D");
      recontree->Branch("y",y,"y[ncamera][15]/D");
      recontree->Branch("ntracks",ntracks,"ntrac ks[ncamera]/I");
      recontree->Branch("edge",edge,"edge[ncamera][15]/O");
      recontree->Branch("skewness",skewness,"skewness[ncamera][15]/D");
      recontree->Branch("date",&date,"date/I");
      recontree->Branch("time",&time,"time/I");
      gROOT->cd();

      TClonesArray* cleanimage = 0;
      bool spark[ncamera];
      double integral[ncamera];
      rawtree->SetBranchAddress("cleanimage",&cleanimage);
      rawtree->SetBranchAddress("spark",spark);
      rawtree->SetBranchAddress("integral",integral);

      MaxCamCluster* clust[15];

      //Add a skimming file
      TString skimfilename = routfilename;
      skimfilename.ReplaceAll(".root","skim.root");
      TFile* skimfile = new TFile(recondir+skimfilename, "RECREATE");
      skimfile->cd();
      TTree* skimtree = new TTree("skim","Skimmed Events");

      int eventnum;
      TClonesArray* skimcleanimage =  new TClonesArray("TH2F",2);
      bool skimspark[ncamera];
      double skimintegral[ncamera];
      DmtpcEvent* skimevent;
      TObjArray* clustarr = new TObjArray(15);
      
      skimtree->Branch("event","DmtpcEvent",&skimevent,128000,0);
      skimtree->Branch("ncamera",&ncam,"ncamera/I");
      skimtree->Branch("cleanimage","TClonesArray",&skimcleanimage,128000,0);
      skimtree->Branch("spark",&skimspark,"spark[ncamera]/O");
      skimtree->Branch("integral",&skimintegral,"integral[ncamera]/D");
      skimtree->Branch("eventnum",&eventnum,"eventnum/I");
      skimtree->Branch("theta",theta,"theta[ncamera][15]/D");
      skimtree->Branch("phi",phi,"phi[ncamera][15]/D");
      skimtree->Branch("E",E,"E[ncamera][15]/D");
      skimtree->Branch("range",range,"range[ncamera][15]/D");
      skimtree->Branch("x",x,"x[ncamera][15]/D");
      skimtree->Branch("y",y,"y[ncamera][15]/D");
      skimtree->Branch("ntracks",ntracks,"ntracks[ncamera]/I");
      skimtree->Branch("edge",edge,"edge[ncamera][15]/O");
      skimtree->Branch("skewness",skewness,"skewness[ncamera][15]/D");
      skimtree->Branch("date",&date,"date/I");
      skimtree->Branch("time",&time,"time/I");
      skimtree->Branch("clusters","TObjArray",&clustarr,128000,0);

      gROOT->cd();




      TCanvas* c  = new TCanvas("c","c",0,0,1000,500);
      c->Divide(2,1);

      int nempty=0;

      //event selection and reconstruction should go here
//      for(int i = 0; i<rawtree->GetEntries(); i++)
      for(int i = 0; i<300; i++) //for testing
      {
	 cout << i << endl;
	 d.getEvent(i);
	 rawtree->GetEvent(i);
	 date = d.event()->timeStamp()->GetDate();
	 time = d.event()->timeStamp()->GetTime();


	 int sumtracks=0;
	 for(int u =0; u<ncamera; u++)
	 {
	    skimspark[u]=spark[u];
	    skimintegral[u]=integral[u];
	    //cloned image
	    TH2F* baseimage = (TH2F*)cleanimage->At(u);
	    //rebin, blur, find cluster
	    if(spark[u]==0)
	    {
	       baseimage->Rebin2D(2,2);
	       baseimage = MaxCamImageTools::blur(baseimage,1,0.8);
	       ntracks[u]=MaxCamImageTools::findClusters(baseimage,
							 clust,18,100000,
							 11,1600);
	       sumtracks+=ntracks[u];
	       
	       TH2F* clustimg = (TH2F*)baseimage;
	       
	       for(int v =0; v<ntracks[u]; v++)
	       {
		  clustarr->Add(clust[v]);
		  
		  c->cd(1);
		  baseimage->DrawCopy("colz");
		  vector<int> bins = clust[v]->getCluster();
		  for(int k=0; k<int(bins.size()); k++)
		  {
		     clustimg->SetBinContent(bins[k],100*(v+1));
		  }
		  c->cd(2);
		  clustimg->Draw("colz");		
		  c->Update();
		  
		  //get reconstruction variables
		  theta[u][v] = 0;
		  phi[u][v] = clust[v]->getPhi2();
		  E[u][v] = clust[v]->getIntegral();
		  range[u][v] = clust[v]->getLength2(phi[u][v],1);
		  clust[v]->getXY(x[u][v],y[u][v]);
		  edge[u][v] = clust[v]->hitsEdge();
		  skewness[u][v] = clust[v]->getSkewness(phi[u][v]);

	       }
	       for(int v=ntracks[u]; v<15; v++)
	       {
		  theta[u][v] = 0;
		  phi[u][v]=0;
		  E[u][v]=0;
		  range[u][v]=0;
		  x[u][v]=0;
		  y[u][v]=0;
		  edge[u][v]=0;
		  skewness[u][v] = 0;
	       }
	    }
	    else
	    {   
	       ntracks[u]=0;
	    }
	 }
	 routfile->cd();
	 recontree->Fill();
	 gROOT->cd();
	 
	 if(sumtracks == 0) nempty++;
	 if(sumtracks>0 || nempty%100 ==0)
	 {
	    eventnum=i;
	    skimevent = d.event();
	    skimcleanimage = cleanimage;
	    skimfile->cd();
	    skimtree->Fill();
	    gROOT->cd();
	 }

      }
     routfile->cd();
     recontree->Write();
     gROOT->cd();
     
     delete recontree;
     routfile->Close();

     skimfile->cd();
     skimtree->Write();
     gROOT->cd();

     delete skimtree;
     skimfile->Close();
     
     
   }

   return 0;
   
}
