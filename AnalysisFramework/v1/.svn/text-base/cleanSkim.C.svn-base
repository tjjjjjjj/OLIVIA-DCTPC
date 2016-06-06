
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


int cleanSkim(TString rawdatafiles="files.txt",TString keyfilename = "keys.txt", TString outdir="./skim/")
{

   TString key = "skim";

   const int nreqmod = 0;
   TString reqmod[nreqmod];
   
   bool pass;
   DmtpcKeys k(keyfilename,rawdatafiles,key,outdir,nreqmod,reqmod,pass);
   if(!pass) return -1;

   TCanvas* c = new TCanvas("c","c",0,0,800,800);
   c->Divide(2,2);
   TCanvas* cc = new TCanvas("cc","cc",850,0,800,400);
   cc->Divide(2,1);


   TFile* routfile;
   TTree* skimtree;
   TH2F* tempimg;

   for(int f=0; f<k.getNFiles(); f++)
   {
      cout << k.getFile(f) << endl;
      
      //Create DMPTC Database and draw out tree
      DmtpcDataset d;
      d.openRootFile(k.getRootDirName()+k.getFile(f));
      TTree* rawtree = k.getBaseTree(f);
      
      int runnum;
      TString runnums = k.getFile(f);
      runnums.ReplaceAll("dmtpc_run","");
      runnums.ReplaceAll(".root","");
      if(!(runnums.IsDigit()))
      {
	 cout << "cannot extract run number" << endl;
	 runnum=0;
      }
      else
	 runnum = runnums.Atoi();

      //Add friends
      k.addFriends(rawtree,f);

      //name outputfile, open up outputfile
      TString routfilename = k.getFile(f);
      routfilename.ReplaceAll(".root",key+".root");
      routfile = new TFile(outdir+routfilename,"CREATE");
      
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
      int ncam = ncamera;

      //set sigma cut values this should be database variable, 
      //find out how to extract for multicamera runs
      double sigmathr;
      if(runnum < 326)
	 sigmathr = 9;
      if(runnum >=326)
	 sigmathr = 20;

      //create skim tree to store reconstruction variables
      skimtree = new TTree(key,"Skimmed Events");

      double theta[ncamera][15], phi[ncamera][15], E[ncamera][15];
      double range[ncamera][15], x[ncamera][15], y[ncamera][15];
      double skewness[ncamera][15];
      bool edge[ncamera][15];
      int ntracks[ncamera];
      int date;
      int time;
      int eventnum;
      TClonesArray* cleanimage = new TClonesArray("TH2F",ncamera);
      bool spark[ncamera];
      double integral[ncamera];      
      DmtpcEvent* skimevent = new DmtpcEvent();
      TObjArray* clustarr = new TObjArray(ncamera);
      TObjArray* cl = new TObjArray(15);
      MaxCamCluster* clust[15];


      skimtree->Branch("event","DmtpcEvent",&skimevent,128000,0);
      skimtree->Branch("ncamera",&ncam,"ncamera/I");
      skimtree->Branch("cleanimage","TClonesArray",&cleanimage,128000,0);
      skimtree->Branch("spark",&spark,"spark[ncamera]/O");
      skimtree->Branch("integral",&integral,"integral[ncamera]/D");
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

      //open up file for bias frame saving
      TString biasoutfilename = routfilename;
      biasoutfilename.ReplaceAll(key+".root","bias.root");

      TFile* biasoutfile = new TFile(outdir+biasoutfilename,"RECREATE");

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
//	 sparkcut[j] =  MaxCamImageTools::countPixelsAboveThreshold(d.getBiasFrame(j+1),countthresh[j])-2000;
	 sparkcut[j]=46000;

	 cout << "stdev: " << MaxCamImageTools::getRMS(d.getBiasFrame(j+1));


	 cout << "\t count: " << countthresh[j] << "\t cut: " << sparkcut[j] << endl;
	 
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

	    if(nabove<sparkcut[u]
	       && MaxCamImageTools::getRMS(tempimg) < sigmathr)
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

      cout << nframes[0] << endl;

      for(int i=0; i<ncamera; i++)
      {
	 if(nframes[i] !=0 && runnum < 500)
	    secondarybias[i]->Scale(1/double(nframes[i]));
	 else
	    secondarybias[i] = (TH2F*)d.getBiasFrame(i+1)->Clone();

	 MaxCamImageTools::killLonePixels2(secondarybias[i],1500);
	 cc->cd(i+1);
	 secondarybias[i]->Draw("colz");

 	 biasframe = secondarybias[i];
 	 biasoutfile->cd();
 	 biastree->Fill();
 	 gROOT->cd();
	 
      }
      cc->Update();

      cout << "All preclean activities done" << endl;
      int nempty=0;

      //cleaning should go here
      for(int i = 0; i<nev; i++)
//      for(int i = 0; i<50; i++) //For testing
      {
	 cout << i << endl;
	 d.getEvent(i);
	 int sumtracks=0;
	 int sparks=0;

	 clustarr->Clear();

	 for(int u=0; u<ncamera; u++)
	 {
	    cout << "\t" << u << ":" << endl;

	    //create cloned image
	    tempimg = (TH2F*)d.event()->ccdData(u);
	    c->cd(1);
	    tempimg->DrawCopy("colz");

	    //check for sparking
	    int nabove = MaxCamImageTools::countPixelsAboveThreshold(tempimg,countthresh[u]);

	    if(nabove < sparkcut[u])
	       spark[u]=0;
	    else
	    {
	       spark[u]=1;
	       cout << nabove << endl;
	    }
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

	    c->cd(2);
	    tempimg->DrawCopy("colz");
	    
	    if(spark[u]==1) sparks++;
	    cl->Clear();
	    TH2F* baseimage = (TH2F*)tempimg->Clone();
	    //rebin, blur, find cluster
	    if(spark[u]==0)
	    {
 	       baseimage->Rebin2D(2,2);
	       baseimage = MaxCamImageTools::blur(baseimage,1,0.8);
	       ntracks[u]=MaxCamImageTools::findClusters(baseimage,
							 clust,3.7,300,
							 5,3600);
// 	       ntracks[u]=MaxCamImageTools::findClusterSE(baseimage,
// 							 clust,18,100000,
// // 							 11,1600);
	       c->cd(3);
	       baseimage->DrawCopy("colz");

	       sumtracks+=ntracks[u];

	       TH2F* clustimg = (TH2F*)tempimg->Clone();
	       for(int v=0; v<ntracks[u]; v++)
	       {
		  clust[v]->changeImage(clustimg);
		  
		  //Correct for tracks in pedestal subtraction
		  if(v==0)
		  {
		     double thirdped = MaxCamImageTools::findPedestalWithTracks(clustimg,clust,ntracks[u]);
		     MaxCamImageTools::subtractPedestal(tempimg,thirdped);
		     c->cd(4);		       
		     tempimg->DrawCopy("colz");
		  }

		  vector<int> pxs = clust[v]->getCluster();
		  for(int q=0; q<int(pxs.size()); q++)
		  {
		     //     clustimg->SetBinContent(pxs[q],100*(v+1));
		  }
	    
		  cl->Add(clust[v]);
		  
		  //get reconstruction variables
		  theta[u][v] = 0;
		  phi[u][v] = clust[v]->getPhi2();
		  E[u][v] = clust[v]->getIntegral();
		  //double xb,yb,xe,ye;
		  //range[u][v] = clust[v]->getLength(xb,yb,xe,ye);
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
	       c->cd(4);
	       clustimg->Draw("colz");

	    }
	    else
	    {   
	       ntracks[u]=0;
	    }



	    new((*cleanimage)[u]) TH2F(*tempimg);
            delete baseimg;
	    c->Update();
	    TString temp;
// 	    if(ntracks[u]>0)
// 	    cin >> temp;

	    clustarr->Add(cl);
	 }
	 
	 if(sumtracks == 0 && sparks==0) nempty++;
	 if(sumtracks>0 || nempty%100 ==0 || sparks>0)
	 {
	    eventnum=i;
	    skimevent = d.event();
	    routfile->cd();
	    skimtree->Fill();
	    gROOT->cd();
	 }

      }

      routfile->cd();
      skimtree->Write();

      delete skimtree;

      routfile->Close();

      biasoutfile->cd();
      biastree->Write();
      
      delete biastree;
      
      biasoutfile->Close();

   }

      
   return 0;
}
