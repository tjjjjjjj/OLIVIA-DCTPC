
#include "../../../MaxCam/MaxCamMC.hh"
#include "../../../MaxCam/MaxCamImageTools.hh"
#include "../../../MaxCam/MaxCamCluster.hh"
#include "../../../MaxCam/DmtpcEvent.hh"
#include "../../../MaxCam/DmtpcDataset.hh"
#include "../../../MaxCam/DmtpcKeys.hh"
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
#include <string>

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
   TCanvas* cc = new TCanvas("cc","cc",850,0,800,800);
   cc->Divide(2,2);


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
      const int eventVersion = d.event()->IsA()->GetClassVersion();


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
      bool spark[ncamera];
      double integral[ncamera];      
      DmtpcEvent* skimevent = new DmtpcEvent();
      TObjArray* clust = new TObjArray(ncamera);
      MaxCamClusterImage* clusti;

      skimtree->Branch("event","DmtpcEvent",&skimevent,128000,0);
      skimtree->Branch("ncamera",&ncam,"ncamera/I");
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
      skimtree->Branch("clusters","TObjArray",&clust,128000,0);
      gROOT->cd();

      //open up file for bias frame saving
      TString biasoutfilename = routfilename;
      biasoutfilename.ReplaceAll(key+".root","bias.root");
      TFile* biasoutfile = new TFile(outdir+biasoutfilename,"RECREATE");
      TH2F* biasframe = d.event()->ccdData(0);
      TTree* biastree = new TTree("bias","Bias Information");
      biastree->Branch("biasframe","TH2F",&biasframe,128000,0);
      gROOT->cd();

      double sigmathr[ncamera];
      int nframes[ncamera];
      for(int u=0; u<ncamera; u++){nframes[u]=0;}
      TH2F* secondarybias[ncamera];
      double bias_mean[ncamera], bias_rms[ncamera];
      double last_mean[ncamera]; // for the spark cut

      // handle the bias frame
      //      define the bias threshold for runs<500 and MC
      for(int u=0; u<ncamera; u++)
	{
	  cout << "Camera " << u << " image cleaning:" << endl;
	  secondarybias[u] = (TH2F*)d.getBiasFrame(u+1)->Clone();
	  cc->cd(u+1);
	  secondarybias[u]->Draw("colz");
	  cc->Update();

	  // the cluster finding threshold is very sensitive to this number!
	  // using 1.5 sigma reproduces cluster finding efficiency from AF v1
	  sigmathr[u]=1.5*MaxCamImageTools::getRMS(secondarybias[u]);
	  cout << "     camera bias MEAN: " << 
	       MaxCamImageTools::getMean(secondarybias[u]) << endl;
	  cout << "     camera bias RMS: " << 
	       MaxCamImageTools::getRMS(secondarybias[u]) << endl;
	}

      // handle the bias frame
      const int nev = rawtree->GetEntries();
      // grab every 100th image to make bias frame; check for sparks
      if ( runnum < 500 ) {
	for(int i=0; i<nev; i=i+100)
	  {
	    cout << "making bias frame: event " << i << endl;
	    d.getEvent(i);
	    for(int u=0; u<ncamera; u++)
	      {
		tempimg = (TH2F*)d.event()->ccdData(u)->Clone("backimg");
		// clean up the bias frame
		MaxCamImageTools::meanRMSNoOutliers(tempimg, bias_mean[u], bias_rms[u]);
		cout << "     bias image mean (no outliers): " << bias_mean[u] << endl;
		cout << "     bias image rms (no outliers): " << bias_rms[u] << endl;
		// make sure the bias is reasonable
		if( (MaxCamImageTools::getRMS(tempimg) < sigmathr[u]) )
		  {
		    // kill 5 sigma outliers in the bias frame, iterate 3 times
		    for (int it=0; it<3; it++) {
		      double mean, rms;
		      mean=MaxCamImageTools::getMean(tempimg);
		      rms=MaxCamImageTools::getRMS(tempimg);
		      if (runnum == 0) { // MC
			  MaxCamImageTools::killPixels(tempimg,mean+10.*rms);
		      } else {
			  MaxCamImageTools::killPixels(tempimg,mean*1.25);
		      }
		      bias_mean[u]=mean;
		      bias_rms[u]=rms;
		    }
		    if(nframes[u]==0) 
		      secondarybias[u]=(TH2F*)tempimg->Clone("average");
		    else
		      secondarybias[u]->Add(tempimg,1);
		    nframes[u]++;
		  }
	      }
	  }
	gROOT->Delete("backimg;*");
      } // end if runnum<500

      for(int u=0; u<ncamera; u++)
	{
	  cout << "number of frames averaged to make camera " << u << " bias: " << nframes[u] << endl;
	  if(nframes[u] !=0 && runnum < 500) {
	    cout << "run number " << runnum << " using averaged bias " << endl;
	    secondarybias[u]->Scale(1/double(nframes[u]));
	  } else {
	    secondarybias[u] = (TH2F*)d.getBiasFrame(u+1)->Clone();
	    cout << "run number " << runnum << " using camera bias " << endl;
	  }
	  // kill 5 sigma outliers in the bias frame, iterate 3 times
	  for (int it=0; it<3; it++) {
	    double mean, rms;
 	    mean=MaxCamImageTools::getMean(secondarybias[u]);
	    rms=MaxCamImageTools::getRMS(secondarybias[u]);
	    if (runnum == 0) { // MC
	      MaxCamImageTools::killPixels(secondarybias[u],mean+10.*rms);
	    } else {
	      MaxCamImageTools::killPixels(secondarybias[u],mean*1.25);
	    }
	    cout << "     iteration " << it << " averaged bias image mean (no outliers): " << mean << endl;
	    cout << "     iteration " << it << " averaged bias image rms (no outliers): " << rms << endl;
	    bias_mean[u]=mean;
	    bias_rms[u]=rms;
	  }
	  cc->cd(u+3);
	  secondarybias[u]->Draw("colz");	  
	  biasframe = secondarybias[u];
	  biasoutfile->cd();
	  biastree->Fill();
	  gROOT->cd();	 
	}
      cc->Update();
      cout << "All preclean activities done" << endl;
      int nempty=0;


      // loop over events
      for(int i = 0; i<nev; i++)
	 //for(int i = 0; i<1000; i++) //For testing
      {
	 cout << i << endl;
	 d.getEvent(i);
	 int sumtracks=0;
	 int sparks=0;
	 
	 clust->Clear();

	 for(int u=0; u<ncamera; u++)
	 {
	    cout << "\t" << u << ":" << endl;
	    // keep track of the lat frame for the spark cut
	    if (i==0) {
	      last_mean[u] = bias_mean[u];
	    }


	    // create cloned image
	    tempimg = (TH2F*)d.event()->ccdData(u)->Clone("tempimg");
	    c->cd(1);
	    tempimg->DrawCopy("colz");

	    double mean, rms;
 	    mean=MaxCamImageTools::getMean(tempimg);
 	    rms=MaxCamImageTools::getRMS(tempimg);
	    double mean_nooutliers, rms_nooutliers;
	    MaxCamImageTools::meanRMSNoOutliers(tempimg, mean_nooutliers, rms_nooutliers);	    
 	    cout << "     image raw mean: " << mean << endl;
	    cout << "     image raw rms: " << rms << endl;
	    cout << "     image n.o. mean: " << mean_nooutliers << endl;
	    cout << "     image n.o. rms: " << rms_nooutliers << endl;
 	    cout << "     bias mean: " << bias_mean[u] << endl;
	    cout << "     bias rms: " << bias_rms[u] << endl;
	    if (runnum > 0) 
	      cout << "     spark test: " << (mean/last_mean[u]) << endl;

	    // check for sparks
	    if( (mean/last_mean[u])<1.01) {
	      spark[u]=0;
	      last_mean[u]=mean;
	    } else {
	      spark[u]=1;
	    }
	    if (runnum == 0) spark[u] = 0; // MC, for now
	    if(spark[u] == 1) cout << "     SPARK!" << endl;

	    // clean up image before cluster finding, same way as bias frame
	    if (runnum==0) { // MC
	      MaxCamImageTools::killLonePixels2(tempimg,mean+10.*rms);
	    } else {
	      MaxCamImageTools::killLonePixels2(tempimg,1.25*mean);
	    }
	    //      subtract off bias frame
	    tempimg->Add(secondarybias[u],-1);
	    //      subtract off remaining pedestal
	    double perpx = tempimg->Integral()/(nbinsx*nbinsy);
	    cout << "     image pedestal correction: " << perpx << endl;
	    if(spark[u]==0)
	       MaxCamImageTools::subtractPedestal(tempimg,perpx);
	    integral[u]=rms; // this variable is written to the output tree
	                     // useful to use rms here since the spark cut
                             // now cuts on rms (instead of integral) 
	    c->cd(2);
	    tempimg->DrawCopy("colz");

	    // save spark images	    
	    if(spark[u]==1) 
	    {
	       sparks++;
	       if(eventVersion==1)
	       {
		  clusti=new MaxCamClusterImage(tempimg,d.event()->timeStamp());
	       }
	       else
		  clusti= new MaxCamClusterImage(tempimg,d.event()->UTCtimeStamp());
	    }

	   

	    // find clusters
	    TH2F* baseimage = (TH2F*)tempimg->Clone("baseimage");
	    
	    if(spark[u]==0)
	      {
		// rebin
		baseimage->Rebin2D(2,2);
		// blur
		baseimage = MaxCamImageTools::blur(baseimage,1,0.8);
		// look for clusters
		if(eventVersion==1){
		clusti = new MaxCamClusterImage(baseimage,d.event()->timeStamp());
		}
		else
		   clusti = new MaxCamClusterImage(baseimage,d.event()->UTCtimeStamp());
		ntracks[u]=MaxCamImageTools::findClustersCI(baseimage,
							    clusti,3.7,300,
							    5,3600);
		c->cd(3);
		baseimage->DrawCopy("colz");
		sumtracks+=ntracks[u];
		
		// get reconstructed quantities
		TH2F* clustimg = (TH2F*)tempimg->Clone("clustimg");
		
		cout << "     unsmeared cluster pixel threshold: " << 3.7*rms_nooutliers << endl;
		clusti->changeImageWithThreshold(clustimg,3.7*rms_nooutliers);	
		for(int v=0; v<ntracks[u]; v++)
		  {
		    vector<int> pxs_smear = clusti->getCluster(v);
		    cout << "     smeared cluster size: " << 
		      int(pxs_smear.size()) << endl;
		    vector<int> pxs_nosmear = clusti->getClusterRed(v);
		    cout << "     unsmeared cluster size: " << 
		      int(pxs_nosmear.size()) << endl;
		    E[u][v] = clusti->getIntegral(v);
		    edge[u][v] = clusti->hitsEdge(v);
		    skewness[u][v] = clusti->getSkewness(v,phi[u][v]);
		    theta[u][v] = 0;
		    phi[u][v] = clusti->getPhi2(v);
		    double xb,yb,xe,ye;
		    range[u][v] = clusti->getLength(v,xb,yb,xe,ye);
		    //		  range[u][v] = clusti->getLength2(v,phi[u][v],1);
		    clusti->getXY(v,x[u][v],y[u][v]);

		    // for debugging cluster finding
		    //		    for(int q=0; q<int(pxs_smear.size()); q++)
		    //		      {
		    //			clustimg->SetBinContent(pxs_smear[q],100*(v+1));
		    //		      }

		  }
		//    initialize the rest of the event array
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
		clustimg->DrawCopy("colz");
	      }
	    else
	      {   
		ntracks[u]=0;
	      }
	    
	    	    
	    c->Update();
	    TString temp;
	    // 	    if(ntracks[u]>0 || spark[u] > 0)
	    //	    if(spark[u] >0)
	    //	     	       cin >> temp;
	    clust->Add(clusti);
	 }
	 
	 // fill tree
	 if(sumtracks == 0 && sparks==0) nempty++;
	 if(sumtracks>0 || nempty%100 ==0 || sparks>0)
	 {
	    eventnum=i;
	    skimevent = d.event();
	    routfile->cd();
	    skimtree->Fill();
	    gROOT->cd();
	 }

	 gROOT->Delete("copy");
	 gROOT->Delete("baseimage");
	 gROOT->Delete("clustimg");
	 gROOT->Delete("tempimg");
	 

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

int main(int argn, char ** argv)
{
   if (argn > 4)
   {
      cout << "./clean_skim [rawdatafiles] [keyfilename] [outdir]" << endl;
      exit(1); 
   } 


   //Weird construction here to always use default arguments for cleanSkim
   if (argn>1) 
   {
     TString datafiles = TString(argv[1]);
     if (argn>2)
     {
       TString keyfiles = TString(argv[2]);
       if (argn>3) 
       {
         TString outdir = TString(argv[3]);
         return cleanSkim(datafiles,keyfiles,outdir);
       }
       return cleanSkim(datafiles,keyfiles);
     }
     return cleanSkim(datafiles);
   }
   return cleanSkim(); 
}

